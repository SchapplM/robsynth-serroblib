% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR5
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
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-24 10:36
% Revision: 5d02717ba55fba3c5445be8d9f6bf09c2cd6665f (2019-10-14)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5PRRRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->4), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->4)
	t35 = pkin(9) + qJ(2);
	t37 = qJD(2) * sin(t35);
	t36 = qJD(2) * cos(t35);
	t1 = [0, -t36, 0, 0, 0; 0, -t37, 0, 0, 0; 0, 0, 0, 0, 0; 0, t37, 0, 0, 0; 0, -t36, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:58
	% EndTime: 2019-10-24 10:36:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (29->10), mult. (36->14), div. (0->0), fcn. (36->4), ass. (0->14)
	t41 = sin(qJ(3));
	t46 = qJD(2) * t41;
	t42 = cos(qJ(3));
	t45 = qJD(2) * t42;
	t44 = qJD(3) * t41;
	t43 = qJD(3) * t42;
	t40 = pkin(9) + qJ(2);
	t39 = cos(t40);
	t38 = sin(t40);
	t37 = t38 * t44 - t39 * t45;
	t36 = t38 * t43 + t39 * t46;
	t35 = t38 * t45 + t39 * t44;
	t34 = t38 * t46 - t39 * t43;
	t1 = [0, t37, t34, 0, 0; 0, -t35, -t36, 0, 0; 0, 0, -t44, 0, 0; 0, t36, t35, 0, 0; 0, t34, t37, 0, 0; 0, 0, -t43, 0, 0; 0, -qJD(2) * t38, 0, 0, 0; 0, qJD(2) * t39, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:57
	% EndTime: 2019-10-24 10:36:57
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (87->15), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t77 = qJ(3) + qJ(4);
	t73 = sin(t77);
	t76 = qJD(3) + qJD(4);
	t81 = t76 * t73;
	t74 = cos(t77);
	t80 = t76 * t74;
	t79 = qJD(2) * t73;
	t78 = qJD(2) * t74;
	t75 = pkin(9) + qJ(2);
	t72 = cos(t75);
	t71 = sin(t75);
	t70 = t71 * t81 - t72 * t78;
	t69 = t71 * t80 + t72 * t79;
	t68 = t71 * t78 + t72 * t81;
	t67 = t71 * t79 - t72 * t80;
	t1 = [0, t70, t67, t67, 0; 0, -t68, -t69, -t69, 0; 0, 0, -t81, -t81, 0; 0, t69, t68, t68, 0; 0, t67, t70, t70, 0; 0, 0, -t80, -t80, 0; 0, -qJD(2) * t71, 0, 0, 0; 0, qJD(2) * t72, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-24 10:36:58
	% EndTime: 2019-10-24 10:36:58
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (177->18), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->16)
	t97 = qJ(3) + qJ(4) + qJ(5);
	t92 = sin(t97);
	t96 = qJD(3) + qJD(4) + qJD(5);
	t102 = t96 * t92;
	t93 = cos(t97);
	t101 = t96 * t93;
	t98 = pkin(9) + qJ(2);
	t94 = sin(t98);
	t100 = qJD(2) * t94;
	t95 = cos(t98);
	t99 = qJD(2) * t95;
	t91 = t94 * t102 - t93 * t99;
	t90 = t94 * t101 + t92 * t99;
	t89 = t93 * t100 + t95 * t102;
	t88 = t92 * t100 - t95 * t101;
	t1 = [0, t91, t88, t88, t88; 0, -t89, -t90, -t90, -t90; 0, 0, -t102, -t102, -t102; 0, t90, t89, t89, t89; 0, t88, t91, t91, t91; 0, 0, -t101, -t101, -t101; 0, -t100, 0, 0, 0; 0, t99, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end