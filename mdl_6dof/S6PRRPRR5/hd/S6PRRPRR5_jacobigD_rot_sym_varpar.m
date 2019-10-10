% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 22:33
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(6)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(11));
	t79 = sin(pkin(11));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t82, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t116 = sin(qJ(2));
	t118 = cos(pkin(6)) * t116;
	t117 = cos(qJ(2));
	t114 = cos(pkin(11));
	t113 = sin(pkin(11));
	t1 = [0, 0, (-t113 * t118 + t114 * t117) * qJD(2), 0, 0, 0; 0, 0, (t113 * t117 + t114 * t118) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t116, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:25
	% EndTime: 2019-10-09 22:33:25
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (12->10), mult. (48->29), div. (0->0), fcn. (50->8), ass. (0->16)
	t133 = sin(pkin(6));
	t136 = sin(qJ(3));
	t146 = t133 * t136;
	t137 = sin(qJ(2));
	t145 = t133 * t137;
	t135 = cos(pkin(6));
	t144 = t135 * t137;
	t139 = cos(qJ(2));
	t143 = t135 * t139;
	t142 = qJD(2) * t136;
	t132 = sin(pkin(11));
	t134 = cos(pkin(11));
	t141 = t132 * t139 + t134 * t144;
	t140 = -t132 * t144 + t134 * t139;
	t138 = cos(qJ(3));
	t1 = [0, 0, t140 * qJD(2), 0, (t132 * t146 + t140 * t138) * qJD(3) + (-t132 * t143 - t134 * t137) * t142, 0; 0, 0, t141 * qJD(2), 0, (-t134 * t146 + t141 * t138) * qJD(3) + (-t132 * t137 + t134 * t143) * t142, 0; 0, 0, qJD(2) * t145, 0, t133 * t139 * t142 + (t135 * t136 + t138 * t145) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:33:26
	% EndTime: 2019-10-09 22:33:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
	t158 = sin(pkin(6));
	t161 = sin(qJ(3));
	t171 = t158 * t161;
	t162 = sin(qJ(2));
	t170 = t158 * t162;
	t160 = cos(pkin(6));
	t169 = t160 * t162;
	t164 = cos(qJ(2));
	t168 = t160 * t164;
	t167 = qJD(2) * t161;
	t157 = sin(pkin(11));
	t159 = cos(pkin(11));
	t166 = t157 * t164 + t159 * t169;
	t165 = -t157 * t169 + t159 * t164;
	t163 = cos(qJ(3));
	t156 = t158 * t164 * t167 + (t160 * t161 + t163 * t170) * qJD(3);
	t155 = (t157 * t171 + t165 * t163) * qJD(3) + (-t157 * t168 - t159 * t162) * t167;
	t154 = (-t159 * t171 + t166 * t163) * qJD(3) + (-t157 * t162 + t159 * t168) * t167;
	t1 = [0, 0, t165 * qJD(2), 0, t155, t155; 0, 0, t166 * qJD(2), 0, t154, t154; 0, 0, qJD(2) * t170, 0, t156, t156;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end