% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RPRRP2
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
%   Siehe auch: S5RPRRP2_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RPRRP2_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRRP2_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_jacobiRD_rot_sym_varpar: pkin has to be [8x1] (double)');
JRD_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:28:31
	% EndTime: 2022-01-23 09:28:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:28:30
	% EndTime: 2022-01-23 09:28:31
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
	% StartTime: 2022-01-23 09:28:31
	% EndTime: 2022-01-23 09:28:31
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
	% StartTime: 2022-01-23 09:28:31
	% EndTime: 2022-01-23 09:28:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t48 = qJ(1) + pkin(8) + qJ(3);
	t49 = qJD(1) + qJD(3);
	t50 = t49 * cos(t48);
	t45 = t49 * sin(t48);
	t1 = [-t50, 0, -t50, 0, 0; -t45, 0, -t45, 0, 0; 0, 0, 0, 0, 0; t45, 0, t45, 0, 0; -t50, 0, -t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:28:31
	% EndTime: 2022-01-23 09:28:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t70 = qJ(1) + pkin(8) + qJ(3);
	t68 = sin(t70);
	t71 = qJD(1) + qJD(3);
	t78 = t71 * t68;
	t72 = sin(qJ(4));
	t77 = t71 * t72;
	t73 = cos(qJ(4));
	t76 = t71 * t73;
	t75 = qJD(4) * t72;
	t74 = qJD(4) * t73;
	t69 = cos(t70);
	t67 = t71 * t69;
	t66 = t68 * t75 - t69 * t76;
	t65 = t68 * t74 + t69 * t77;
	t64 = t68 * t76 + t69 * t75;
	t63 = t68 * t77 - t69 * t74;
	t1 = [t66, 0, t66, t63, 0; -t64, 0, -t64, -t65, 0; 0, 0, 0, -t75, 0; t65, 0, t65, t64, 0; t63, 0, t63, t66, 0; 0, 0, 0, -t74, 0; -t78, 0, -t78, 0, 0; t67, 0, t67, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2022-01-23 09:28:31
	% EndTime: 2022-01-23 09:28:31
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (88->13), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t81 = qJ(1) + pkin(8) + qJ(3);
	t79 = sin(t81);
	t82 = qJD(1) + qJD(3);
	t89 = t82 * t79;
	t83 = sin(qJ(4));
	t88 = t82 * t83;
	t84 = cos(qJ(4));
	t87 = t82 * t84;
	t86 = qJD(4) * t83;
	t85 = qJD(4) * t84;
	t80 = cos(t81);
	t78 = t82 * t80;
	t77 = t79 * t86 - t80 * t87;
	t76 = t79 * t85 + t80 * t88;
	t75 = t79 * t87 + t80 * t86;
	t74 = t79 * t88 - t80 * t85;
	t1 = [t77, 0, t77, t74, 0; -t75, 0, -t75, -t76, 0; 0, 0, 0, -t86, 0; t76, 0, t76, t75, 0; t74, 0, t74, t77, 0; 0, 0, 0, -t85, 0; -t89, 0, -t89, 0, 0; t78, 0, t78, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
end