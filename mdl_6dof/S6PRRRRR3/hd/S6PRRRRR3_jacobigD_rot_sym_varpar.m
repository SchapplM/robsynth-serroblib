% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR3
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 23:17
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR3_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR3_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:28
	% EndTime: 2019-10-09 23:17:28
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(12));
	t79 = sin(pkin(12));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t125 = sin(pkin(6));
	t128 = sin(qJ(3));
	t138 = t125 * t128;
	t129 = sin(qJ(2));
	t137 = t125 * t129;
	t127 = cos(pkin(6));
	t136 = t127 * t129;
	t131 = cos(qJ(2));
	t135 = t127 * t131;
	t134 = qJD(2) * t128;
	t124 = sin(pkin(12));
	t126 = cos(pkin(12));
	t133 = t124 * t131 + t126 * t136;
	t132 = -t124 * t136 + t126 * t131;
	t130 = cos(qJ(3));
	t1 = [0, 0, t132 * qJD(2), (t124 * t138 + t132 * t130) * qJD(3) + (-t124 * t135 - t126 * t129) * t134, 0, 0; 0, 0, t133 * qJD(2), (-t126 * t138 + t133 * t130) * qJD(3) + (-t124 * t129 + t126 * t135) * t134, 0, 0; 0, 0, qJD(2) * t137, t125 * t131 * t134 + (t127 * t128 + t130 * t137) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
	t154 = sin(pkin(6));
	t157 = sin(qJ(3));
	t167 = t154 * t157;
	t158 = sin(qJ(2));
	t166 = t154 * t158;
	t156 = cos(pkin(6));
	t165 = t156 * t158;
	t160 = cos(qJ(2));
	t164 = t156 * t160;
	t163 = qJD(2) * t157;
	t153 = sin(pkin(12));
	t155 = cos(pkin(12));
	t162 = t153 * t160 + t155 * t165;
	t161 = -t153 * t165 + t155 * t160;
	t159 = cos(qJ(3));
	t152 = t154 * t160 * t163 + (t156 * t157 + t159 * t166) * qJD(3);
	t151 = (t153 * t167 + t161 * t159) * qJD(3) + (-t153 * t164 - t155 * t158) * t163;
	t150 = (-t155 * t167 + t162 * t159) * qJD(3) + (-t153 * t158 + t155 * t164) * t163;
	t1 = [0, 0, t161 * qJD(2), t151, t151, 0; 0, 0, t162 * qJD(2), t150, t150, 0; 0, 0, qJD(2) * t166, t152, t152, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:17:29
	% EndTime: 2019-10-09 23:17:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (32->10), mult. (120->29), div. (0->0), fcn. (126->8), ass. (0->19)
	t159 = sin(pkin(6));
	t162 = sin(qJ(3));
	t172 = t159 * t162;
	t163 = sin(qJ(2));
	t171 = t159 * t163;
	t161 = cos(pkin(6));
	t170 = t161 * t163;
	t165 = cos(qJ(2));
	t169 = t161 * t165;
	t168 = qJD(2) * t162;
	t158 = sin(pkin(12));
	t160 = cos(pkin(12));
	t167 = t158 * t165 + t160 * t170;
	t166 = -t158 * t170 + t160 * t165;
	t164 = cos(qJ(3));
	t157 = t159 * t165 * t168 + (t161 * t162 + t164 * t171) * qJD(3);
	t156 = (t158 * t172 + t166 * t164) * qJD(3) + (-t158 * t169 - t160 * t163) * t168;
	t155 = (-t160 * t172 + t167 * t164) * qJD(3) + (-t158 * t163 + t160 * t169) * t168;
	t1 = [0, 0, t166 * qJD(2), t156, t156, t156; 0, 0, t167 * qJD(2), t155, t155, t155; 0, 0, qJD(2) * t171, t157, t157, t157;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end