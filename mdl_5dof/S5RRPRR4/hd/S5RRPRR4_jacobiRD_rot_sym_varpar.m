% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRPRR4
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
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRPRR4_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRR4_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
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
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
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
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (30->8), mult. (8->2), div. (0->0), fcn. (8->2), ass. (0->5)
	t29 = qJ(1) + qJ(2) + pkin(9);
	t30 = qJD(1) + qJD(2);
	t31 = t30 * sin(t29);
	t26 = t30 * cos(t29);
	t1 = [0, 0, 0, 0, 0; -t31, -t31, 0, 0, 0; t26, t26, 0, 0, 0; 0, 0, 0, 0, 0; -t26, -t26, 0, 0, 0; -t31, -t31, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (86->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t123 = qJD(1) + qJD(2);
	t124 = sin(qJ(4));
	t129 = t123 * t124;
	t125 = cos(qJ(4));
	t128 = t123 * t125;
	t127 = qJD(4) * t124;
	t126 = qJD(4) * t125;
	t122 = qJ(1) + qJ(2) + pkin(9);
	t121 = cos(t122);
	t120 = sin(t122);
	t119 = t123 * t121;
	t118 = t123 * t120;
	t117 = -t120 * t127 + t121 * t128;
	t116 = -t120 * t126 - t121 * t129;
	t115 = -t120 * t128 - t121 * t127;
	t114 = t120 * t129 - t121 * t126;
	t1 = [0, 0, 0, -t127, 0; t115, t115, 0, t116, 0; t117, t117, 0, -t114, 0; 0, 0, 0, -t126, 0; t114, t114, 0, -t117, 0; t116, t116, 0, t115, 0; 0, 0, 0, 0, 0; t119, t119, 0, 0, 0; t118, t118, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:02:38
	% EndTime: 2020-01-03 12:02:38
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (168->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t160 = qJ(4) + qJ(5);
	t156 = sin(t160);
	t158 = qJD(4) + qJD(5);
	t162 = t158 * t156;
	t157 = cos(t160);
	t161 = t158 * t157;
	t155 = qJ(1) + qJ(2) + pkin(9);
	t153 = sin(t155);
	t159 = qJD(1) + qJD(2);
	t151 = t159 * t153;
	t154 = cos(t155);
	t152 = t159 * t154;
	t148 = t157 * t152 - t153 * t162;
	t147 = -t156 * t152 - t153 * t161;
	t146 = -t157 * t151 - t154 * t162;
	t145 = t156 * t151 - t154 * t161;
	t1 = [0, 0, 0, -t162, -t162; t146, t146, 0, t147, t147; t148, t148, 0, -t145, -t145; 0, 0, 0, -t161, -t161; t145, t145, 0, -t148, -t148; t147, t147, 0, t146, t146; 0, 0, 0, 0, 0; t152, t152, 0, 0, 0; t151, t151, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end