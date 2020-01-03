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
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (3->3), mult. (4->2), div. (0->0), fcn. (4->2), ass. (0->3)
	t10 = qJD(1) * sin(qJ(1));
	t9 = qJD(1) * cos(qJ(1));
	t1 = [0, 0, 0, 0, 0; -t10, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t10, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiRD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (22->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t23 = qJD(1) + qJD(2);
	t24 = qJ(1) + qJ(2);
	t25 = t23 * sin(t24);
	t20 = t23 * cos(t24);
	t1 = [0, 0, 0, 0, 0; -t25, -t25, 0, 0, 0; t20, t20, 0, 0, 0; 0, 0, 0, 0, 0; -t20, -t20, 0, 0, 0; -t25, -t25, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiRD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (28->6), mult. (20->8), div. (0->0), fcn. (20->4), ass. (0->13)
	t83 = qJD(1) + qJD(2);
	t90 = t83 * sin(pkin(9));
	t89 = t83 * cos(pkin(9));
	t84 = qJ(1) + qJ(2);
	t81 = sin(t84);
	t88 = t81 * t89;
	t82 = cos(t84);
	t87 = t82 * t90;
	t80 = t83 * t82;
	t79 = t83 * t81;
	t78 = t82 * t89;
	t77 = t81 * t90;
	t1 = [0, 0, 0, 0, 0; -t88, -t88, 0, 0, 0; t78, t78, 0, 0, 0; 0, 0, 0, 0, 0; t77, t77, 0, 0, 0; -t87, -t87, 0, 0, 0; 0, 0, 0, 0, 0; t80, t80, 0, 0, 0; t79, t79, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (84->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t123 = qJ(1) + qJ(2);
	t119 = sin(t123);
	t122 = qJD(1) + qJD(2);
	t115 = t122 * t119;
	t120 = cos(t123);
	t116 = t122 * t120;
	t125 = qJD(4) * t119;
	t124 = qJD(4) * t120;
	t121 = pkin(9) + qJ(4);
	t118 = cos(t121);
	t117 = sin(t121);
	t114 = t118 * t116 - t117 * t125;
	t113 = -t117 * t116 - t118 * t125;
	t112 = -t118 * t115 - t117 * t124;
	t111 = t117 * t115 - t118 * t124;
	t1 = [0, 0, 0, -qJD(4) * t117, 0; t112, t112, 0, t113, 0; t114, t114, 0, -t111, 0; 0, 0, 0, -qJD(4) * t118, 0; t111, t111, 0, -t114, 0; t113, t113, 0, t112, 0; 0, 0, 0, 0, 0; t116, t116, 0, 0, 0; t115, t115, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:04:36
	% EndTime: 2020-01-03 12:04:37
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (168->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t150 = pkin(9) + qJ(4) + qJ(5);
	t148 = sin(t150);
	t153 = qJD(4) + qJD(5);
	t157 = t153 * t148;
	t149 = cos(t150);
	t156 = t153 * t149;
	t155 = qJ(1) + qJ(2);
	t151 = sin(t155);
	t154 = qJD(1) + qJD(2);
	t146 = t154 * t151;
	t152 = cos(t155);
	t147 = t154 * t152;
	t143 = t149 * t147 - t151 * t157;
	t142 = -t148 * t147 - t151 * t156;
	t141 = -t149 * t146 - t152 * t157;
	t140 = t148 * t146 - t152 * t156;
	t1 = [0, 0, 0, -t157, -t157; t141, t141, 0, t142, t142; t143, t143, 0, -t140, -t140; 0, 0, 0, -t156, -t156; t140, t140, 0, -t143, -t143; t142, t142, 0, t141, t141; 0, 0, 0, 0, 0; t147, t147, 0, 0, 0; t146, t146, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end