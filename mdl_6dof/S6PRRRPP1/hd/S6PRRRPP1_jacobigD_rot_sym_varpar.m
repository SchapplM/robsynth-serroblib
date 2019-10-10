% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRPP1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPP1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRPP1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(10));
	t79 = sin(pkin(10));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
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
	t124 = sin(pkin(10));
	t126 = cos(pkin(10));
	t133 = t124 * t131 + t126 * t136;
	t132 = -t124 * t136 + t126 * t131;
	t130 = cos(qJ(3));
	t1 = [0, 0, t132 * qJD(2), (t124 * t138 + t132 * t130) * qJD(3) + (-t124 * t135 - t126 * t129) * t134, 0, 0; 0, 0, t133 * qJD(2), (-t126 * t138 + t133 * t130) * qJD(3) + (-t124 * t129 + t126 * t135) * t134, 0, 0; 0, 0, qJD(2) * t137, t125 * t131 * t134 + (t127 * t128 + t130 * t137) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t134 = sin(pkin(6));
	t137 = sin(qJ(3));
	t147 = t134 * t137;
	t138 = sin(qJ(2));
	t146 = t134 * t138;
	t136 = cos(pkin(6));
	t145 = t136 * t138;
	t140 = cos(qJ(2));
	t144 = t136 * t140;
	t143 = qJD(2) * t137;
	t133 = sin(pkin(10));
	t135 = cos(pkin(10));
	t142 = t133 * t140 + t135 * t145;
	t141 = -t133 * t145 + t135 * t140;
	t139 = cos(qJ(3));
	t1 = [0, 0, t141 * qJD(2), (t133 * t147 + t141 * t139) * qJD(3) + (-t133 * t144 - t135 * t138) * t143, 0, 0; 0, 0, t142 * qJD(2), (-t135 * t147 + t142 * t139) * qJD(3) + (-t133 * t138 + t135 * t144) * t143, 0, 0; 0, 0, qJD(2) * t146, t134 * t140 * t143 + (t136 * t137 + t139 * t146) * qJD(3), 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:41:08
	% EndTime: 2019-10-09 22:41:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t163 = sin(pkin(6));
	t166 = sin(qJ(3));
	t176 = t163 * t166;
	t167 = sin(qJ(2));
	t175 = t163 * t167;
	t165 = cos(pkin(6));
	t174 = t165 * t167;
	t169 = cos(qJ(2));
	t173 = t165 * t169;
	t172 = qJD(2) * t166;
	t162 = sin(pkin(10));
	t164 = cos(pkin(10));
	t171 = t162 * t169 + t164 * t174;
	t170 = -t162 * t174 + t164 * t169;
	t168 = cos(qJ(3));
	t1 = [0, 0, t170 * qJD(2), (t162 * t176 + t170 * t168) * qJD(3) + (-t162 * t173 - t164 * t167) * t172, 0, 0; 0, 0, t171 * qJD(2), (-t164 * t176 + t171 * t168) * qJD(3) + (-t162 * t167 + t164 * t173) * t172, 0, 0; 0, 0, qJD(2) * t175, t163 * t169 * t172 + (t165 * t166 + t168 * t175) * qJD(3), 0, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end