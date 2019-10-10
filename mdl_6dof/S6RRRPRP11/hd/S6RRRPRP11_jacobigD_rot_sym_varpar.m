% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:51
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:31
	% EndTime: 2019-10-10 11:51:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
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
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t117 = sin(qJ(2));
	t118 = sin(qJ(1));
	t125 = t117 * t118;
	t120 = cos(qJ(1));
	t124 = t117 * t120;
	t119 = cos(qJ(2));
	t123 = t118 * t119;
	t122 = t119 * t120;
	t115 = sin(pkin(6));
	t121 = qJD(1) * t115;
	t116 = cos(pkin(6));
	t1 = [0, t120 * t121, (-t116 * t125 + t122) * qJD(2) + (t116 * t122 - t125) * qJD(1), 0, 0, 0; 0, t118 * t121, (t116 * t124 + t123) * qJD(2) + (t116 * t123 + t124) * qJD(1), 0, 0, 0; 0, 0, t115 * qJD(2) * t117, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t151 = sin(pkin(6));
	t154 = sin(qJ(2));
	t171 = t151 * t154;
	t155 = sin(qJ(1));
	t170 = t151 * t155;
	t158 = cos(qJ(1));
	t169 = t151 * t158;
	t168 = t154 * t155;
	t167 = t154 * t158;
	t157 = cos(qJ(2));
	t166 = t155 * t157;
	t165 = t158 * t157;
	t164 = qJD(1) * t151;
	t156 = cos(qJ(3));
	t163 = qJD(2) * t156;
	t152 = cos(pkin(6));
	t162 = t152 * t165 - t168;
	t161 = t152 * t166 + t167;
	t160 = t152 * t167 + t166;
	t159 = -t152 * t168 + t165;
	t153 = sin(qJ(3));
	t1 = [0, t158 * t164, t162 * qJD(1) + t159 * qJD(2), 0, (-t159 * t153 + t156 * t170) * qJD(3) - t161 * t163 + (t153 * t169 - t160 * t156) * qJD(1), 0; 0, t155 * t164, t161 * qJD(1) + t160 * qJD(2), 0, (-t160 * t153 - t156 * t169) * qJD(3) + t162 * t163 + (t153 * t170 + t159 * t156) * qJD(1), 0; 0, 0, qJD(2) * t171, 0, t151 * t157 * t163 + (t152 * t156 - t153 * t171) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:51:32
	% EndTime: 2019-10-10 11:51:32
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t158 = sin(pkin(6));
	t161 = sin(qJ(2));
	t178 = t158 * t161;
	t162 = sin(qJ(1));
	t177 = t158 * t162;
	t165 = cos(qJ(1));
	t176 = t158 * t165;
	t175 = t161 * t162;
	t174 = t161 * t165;
	t164 = cos(qJ(2));
	t173 = t162 * t164;
	t172 = t165 * t164;
	t171 = qJD(1) * t158;
	t163 = cos(qJ(3));
	t170 = qJD(2) * t163;
	t159 = cos(pkin(6));
	t169 = t159 * t172 - t175;
	t168 = t159 * t173 + t174;
	t167 = t159 * t174 + t173;
	t166 = -t159 * t175 + t172;
	t160 = sin(qJ(3));
	t1 = [0, t165 * t171, t169 * qJD(1) + t166 * qJD(2), 0, (-t166 * t160 + t163 * t177) * qJD(3) - t168 * t170 + (t160 * t176 - t167 * t163) * qJD(1), 0; 0, t162 * t171, t168 * qJD(1) + t167 * qJD(2), 0, (-t167 * t160 - t163 * t176) * qJD(3) + t169 * t170 + (t160 * t177 + t166 * t163) * qJD(1), 0; 0, 0, qJD(2) * t178, 0, t158 * t164 * t170 + (t159 * t163 - t160 * t178) * qJD(3), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end