% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR5_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR5_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:48
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (1->1), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = qJD(1) * cos(qJ(1));
	t7 = qJD(1) * sin(qJ(1));
	t1 = [0, 0, 0, 0, 0; t7, 0, 0, 0, 0; -t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; t9, 0, 0, 0, 0; t7, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (18->4), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t26 = qJ(1) + qJ(2);
	t25 = qJD(1) + qJD(2);
	t23 = t25 * cos(t26);
	t22 = t25 * sin(t26);
	t1 = [0, 0, 0, 0, 0; t22, t22, 0, 0, 0; -t23, -t23, 0, 0, 0; 0, 0, 0, 0, 0; t23, t23, 0, 0, 0; t22, t22, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (32->10), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t75 = qJ(1) + qJ(2);
	t72 = sin(t75);
	t74 = qJD(1) + qJD(2);
	t83 = t74 * t72;
	t73 = cos(t75);
	t82 = t74 * t73;
	t81 = t74 * sin(pkin(9));
	t80 = t74 * cos(pkin(9));
	t79 = t72 * t81;
	t78 = t73 * t80;
	t71 = t73 * t81;
	t70 = t72 * t80;
	t1 = [0, 0, 0, 0, 0; t70, t70, 0, 0, 0; -t78, -t78, 0, 0, 0; 0, 0, 0, 0, 0; -t79, -t79, 0, 0, 0; t71, t71, 0, 0, 0; 0, 0, 0, 0, 0; -t82, -t82, 0, 0, 0; -t83, -t83, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (88->17), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t131 = qJ(1) + qJ(2);
	t127 = sin(t131);
	t130 = qJD(1) + qJD(2);
	t135 = t130 * t127;
	t128 = cos(t131);
	t134 = t130 * t128;
	t133 = qJD(4) * t127;
	t132 = qJD(4) * t128;
	t129 = pkin(9) + qJ(4);
	t126 = cos(t129);
	t125 = sin(t129);
	t122 = -t125 * t133 + t126 * t134;
	t121 = t125 * t134 + t126 * t133;
	t120 = t125 * t132 + t126 * t135;
	t119 = t125 * t135 - t126 * t132;
	t1 = [0, 0, 0, -qJD(4) * t125, 0; t120, t120, 0, t121, 0; -t122, -t122, 0, t119, 0; 0, 0, 0, -qJD(4) * t126, 0; -t119, -t119, 0, t122, 0; t121, t121, 0, t120, 0; 0, 0, 0, 0, 0; -t134, -t134, 0, 0, 0; -t135, -t135, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 18:34:49
	% EndTime: 2019-12-05 18:34:49
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (172->20), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t161 = pkin(9) + qJ(4) + qJ(5);
	t159 = sin(t161);
	t164 = qJD(4) + qJD(5);
	t170 = t164 * t159;
	t160 = cos(t161);
	t169 = t164 * t160;
	t166 = qJ(1) + qJ(2);
	t162 = sin(t166);
	t165 = qJD(1) + qJD(2);
	t168 = t165 * t162;
	t163 = cos(t166);
	t167 = t165 * t163;
	t156 = t160 * t167 - t162 * t170;
	t155 = t159 * t167 + t162 * t169;
	t154 = t160 * t168 + t163 * t170;
	t153 = t159 * t168 - t163 * t169;
	t1 = [0, 0, 0, -t170, -t170; t154, t154, 0, t155, t155; -t156, -t156, 0, t153, t153; 0, 0, 0, -t169, -t169; -t153, -t153, 0, t156, t156; t155, t155, 0, t154, t154; 0, 0, 0, 0, 0; -t167, -t167, 0, 0, 0; -t168, -t168, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end