% Zeitableitung der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% JRD_rot [9x5]
%   Zeitableitung der Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JRD_rot = S5RRRPR3_jacobiRD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR3_jacobiRD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_jacobiRD_rot_sym_varpar: pkin has to be [9x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiRD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiRD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
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
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
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
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:13
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (58->10), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->17)
	t118 = qJD(1) + qJD(2);
	t120 = sin(qJ(3));
	t125 = t118 * t120;
	t121 = cos(qJ(3));
	t124 = t118 * t121;
	t123 = qJD(3) * t120;
	t122 = qJD(3) * t121;
	t119 = qJ(1) + qJ(2);
	t117 = cos(t119);
	t116 = sin(t119);
	t115 = t118 * t117;
	t114 = t118 * t116;
	t113 = -t116 * t123 + t117 * t124;
	t112 = -t116 * t122 - t117 * t125;
	t111 = -t116 * t124 - t117 * t123;
	t110 = t116 * t125 - t117 * t122;
	t1 = [0, 0, -t123, 0, 0; t111, t111, t112, 0, 0; t113, t113, -t110, 0, 0; 0, 0, -t122, 0, 0; t110, t110, -t113, 0, 0; t112, t112, t111, 0, 0; 0, 0, 0, 0, 0; t115, t115, 0, 0, 0; t114, t114, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiRD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:14
	% EndTime: 2020-01-03 12:10:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (84->11), mult. (54->14), div. (0->0), fcn. (54->4), ass. (0->16)
	t135 = qJ(1) + qJ(2);
	t131 = sin(t135);
	t133 = qJD(1) + qJD(2);
	t127 = t133 * t131;
	t132 = cos(t135);
	t128 = t133 * t132;
	t137 = qJD(3) * t131;
	t136 = qJD(3) * t132;
	t134 = qJ(3) + pkin(9);
	t130 = cos(t134);
	t129 = sin(t134);
	t126 = t130 * t128 - t129 * t137;
	t125 = -t129 * t128 - t130 * t137;
	t124 = -t130 * t127 - t129 * t136;
	t123 = t129 * t127 - t130 * t136;
	t1 = [0, 0, -qJD(3) * t129, 0, 0; t124, t124, t125, 0, 0; t126, t126, -t123, 0, 0; 0, 0, -qJD(3) * t130, 0, 0; t123, t123, -t126, 0, 0; t125, t125, t124, 0, 0; 0, 0, 0, 0, 0; t128, t128, 0, 0, 0; t127, t127, 0, 0, 0;];
	JRD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiRD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2020-01-03 12:10:13
	% EndTime: 2020-01-03 12:10:14
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (168->16), mult. (72->12), div. (0->0), fcn. (72->4), ass. (0->17)
	t158 = qJ(3) + pkin(9) + qJ(5);
	t156 = sin(t158);
	t161 = qJD(3) + qJD(5);
	t165 = t161 * t156;
	t157 = cos(t158);
	t164 = t161 * t157;
	t163 = qJ(1) + qJ(2);
	t159 = sin(t163);
	t162 = qJD(1) + qJD(2);
	t154 = t162 * t159;
	t160 = cos(t163);
	t155 = t162 * t160;
	t151 = t157 * t155 - t159 * t165;
	t150 = -t156 * t155 - t159 * t164;
	t149 = -t157 * t154 - t160 * t165;
	t148 = t156 * t154 - t160 * t164;
	t1 = [0, 0, -t165, 0, -t165; t149, t149, t150, 0, t150; t151, t151, -t148, 0, -t148; 0, 0, -t164, 0, -t164; t148, t148, -t151, 0, -t151; t150, t150, t149, 0, t149; 0, 0, 0, 0, 0; t155, t155, 0, 0, 0; t154, t154, 0, 0, 0;];
	JRD_rot = t1;
else
	JRD_rot=NaN(9,5);
end