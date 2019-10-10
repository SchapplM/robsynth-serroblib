% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:29
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPPR8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPPR8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(6));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(6));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t122 = sin(qJ(2));
	t123 = sin(qJ(1));
	t130 = t122 * t123;
	t125 = cos(qJ(1));
	t129 = t122 * t125;
	t124 = cos(qJ(2));
	t128 = t123 * t124;
	t127 = t124 * t125;
	t120 = sin(pkin(6));
	t126 = qJD(1) * t120;
	t121 = cos(pkin(6));
	t1 = [0, t125 * t126, (-t121 * t130 + t127) * qJD(2) + (t121 * t127 - t130) * qJD(1), 0, 0, 0; 0, t123 * t126, (t121 * t129 + t128) * qJD(2) + (t121 * t128 + t129) * qJD(1), 0, 0, 0; 0, 0, t120 * qJD(2) * t122, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:33
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t114 = sin(qJ(2));
	t115 = sin(qJ(1));
	t122 = t114 * t115;
	t117 = cos(qJ(1));
	t121 = t114 * t117;
	t116 = cos(qJ(2));
	t120 = t115 * t116;
	t119 = t116 * t117;
	t112 = sin(pkin(6));
	t118 = qJD(1) * t112;
	t113 = cos(pkin(6));
	t1 = [0, t117 * t118, (-t113 * t122 + t119) * qJD(2) + (t113 * t119 - t122) * qJD(1), 0, 0, 0; 0, t115 * t118, (t113 * t121 + t120) * qJD(2) + (t113 * t120 + t121) * qJD(1), 0, 0, 0; 0, 0, t112 * qJD(2) * t114, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:29:33
	% EndTime: 2019-10-10 11:29:34
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t155 = sin(pkin(6));
	t158 = sin(qJ(2));
	t175 = t155 * t158;
	t159 = sin(qJ(1));
	t174 = t155 * t159;
	t162 = cos(qJ(1));
	t173 = t155 * t162;
	t172 = t158 * t159;
	t171 = t158 * t162;
	t161 = cos(qJ(2));
	t170 = t159 * t161;
	t169 = t162 * t161;
	t168 = qJD(1) * t155;
	t160 = cos(qJ(3));
	t167 = qJD(2) * t160;
	t156 = cos(pkin(6));
	t166 = t156 * t169 - t172;
	t165 = t156 * t170 + t171;
	t164 = t156 * t171 + t170;
	t163 = -t156 * t172 + t169;
	t157 = sin(qJ(3));
	t1 = [0, t162 * t168, qJD(1) * t166 + qJD(2) * t163, 0, 0, (-t157 * t163 + t160 * t174) * qJD(3) - t165 * t167 + (t157 * t173 - t160 * t164) * qJD(1); 0, t159 * t168, qJD(1) * t165 + qJD(2) * t164, 0, 0, (-t157 * t164 - t160 * t173) * qJD(3) + t166 * t167 + (t157 * t174 + t160 * t163) * qJD(1); 0, 0, qJD(2) * t175, 0, 0, t155 * t161 * t167 + (t156 * t160 - t157 * t175) * qJD(3);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end