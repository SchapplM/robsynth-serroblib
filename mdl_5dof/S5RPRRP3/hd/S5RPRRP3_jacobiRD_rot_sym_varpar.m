% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% Zeitableitung: Die Gradientenmatrix wird nochmal nach der Zeit abgeleitet.
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RPRRP3_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t31 = qJD(1) * sin(qJ(1));
	t30 = qJD(1) * cos(qJ(1));
	t1 = [-t30, 0, 0, 0, 0; -t31, 0, 0, 0, 0; 0, 0, 0, 0, 0; t31, 0, 0, 0, 0; -t30, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = qJ(1) + pkin(8);
	t37 = qJD(1) * sin(t35);
	t36 = qJD(1) * cos(t35);
	t1 = [-t36, 0, 0, 0, 0; -t37, 0, 0, 0, 0; 0, 0, 0, 0, 0; t37, 0, 0, 0, 0; -t36, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(1) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(1) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = qJ(1) + pkin(8);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [t37, 0, t34, 0, 0; -t35, 0, -t36, 0, 0; 0, 0, -t44, 0, 0; t36, 0, t35, 0, 0; t34, 0, t37, 0, 0; 0, 0, -t43, 0, 0; -qJD(1) * t38, 0, 0, 0, 0; qJD(1) * t39, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (87->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t77 = qJ(3) + qJ(4);
	t73 = sin(t77);
	t75 = qJD(3) + qJD(4);
	t81 = t75 * t73;
	t74 = cos(t77);
	t80 = t75 * t74;
	t79 = qJD(1) * t73;
	t78 = qJD(1) * t74;
	t76 = qJ(1) + pkin(8);
	t72 = cos(t76);
	t71 = sin(t76);
	t70 = t71 * t81 - t72 * t78;
	t69 = t71 * t80 + t72 * t79;
	t68 = t71 * t78 + t72 * t81;
	t67 = t71 * t79 - t72 * t80;
	t1 = [t70, 0, t67, t67, 0; -t68, 0, -t69, -t69, 0; 0, 0, -t81, -t81, 0; t69, 0, t68, t68, 0; t67, 0, t70, t70, 0; 0, 0, -t80, -t80, 0; -qJD(1) * t71, 0, 0, 0, 0; qJD(1) * t72, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:30:44
	% EndTime: 2022-01-23 09:30:44
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (87->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t86 = qJ(3) + qJ(4);
	t82 = sin(t86);
	t84 = qJD(3) + qJD(4);
	t90 = t84 * t82;
	t83 = cos(t86);
	t89 = t84 * t83;
	t88 = qJD(1) * t82;
	t87 = qJD(1) * t83;
	t85 = qJ(1) + pkin(8);
	t81 = cos(t85);
	t80 = sin(t85);
	t79 = t80 * t90 - t81 * t87;
	t78 = t80 * t89 + t81 * t88;
	t77 = t80 * t87 + t81 * t90;
	t76 = t80 * t88 - t81 * t89;
	t1 = [t79, 0, t76, t76, 0; -t77, 0, -t78, -t78, 0; 0, 0, -t90, -t90, 0; t78, 0, t77, t77, 0; t76, 0, t79, t79, 0; 0, 0, -t89, -t89, 0; -qJD(1) * t80, 0, 0, 0, 0; qJD(1) * t81, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end