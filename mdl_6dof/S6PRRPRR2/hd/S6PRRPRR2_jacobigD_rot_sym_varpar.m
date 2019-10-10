% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR2
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
% Datum: 2019-10-09 22:27
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR2_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR2_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:28
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.03s
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
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t90 = sin(qJ(2));
	t92 = cos(pkin(6)) * t90;
	t91 = cos(qJ(2));
	t88 = cos(pkin(11));
	t87 = sin(pkin(11));
	t1 = [0, 0, (-t87 * t92 + t88 * t91) * qJD(2), 0, 0, 0; 0, 0, (t87 * t91 + t88 * t92) * qJD(2), 0, 0, 0; 0, 0, sin(pkin(6)) * qJD(2) * t90, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (21->11), mult. (48->30), div. (0->0), fcn. (50->8), ass. (0->17)
	t139 = qJ(3) + pkin(12);
	t137 = sin(t139);
	t141 = sin(pkin(6));
	t152 = t141 * t137;
	t143 = cos(pkin(6));
	t144 = sin(qJ(2));
	t151 = t143 * t144;
	t145 = cos(qJ(2));
	t150 = t143 * t145;
	t149 = qJD(2) * t137;
	t148 = qJD(2) * t141;
	t140 = sin(pkin(11));
	t142 = cos(pkin(11));
	t147 = t140 * t145 + t142 * t151;
	t146 = -t140 * t151 + t142 * t145;
	t138 = cos(t139);
	t1 = [0, 0, t146 * qJD(2), 0, (t146 * t138 + t140 * t152) * qJD(3) + (-t140 * t150 - t142 * t144) * t149, 0; 0, 0, t147 * qJD(2), 0, (t147 * t138 - t142 * t152) * qJD(3) + (-t140 * t144 + t142 * t150) * t149, 0; 0, 0, t144 * t148, 0, t145 * t137 * t148 + (t138 * t141 * t144 + t137 * t143) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:27:29
	% EndTime: 2019-10-09 22:27:29
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (40->11), mult. (84->30), div. (0->0), fcn. (88->8), ass. (0->20)
	t167 = qJ(3) + pkin(12);
	t165 = sin(t167);
	t169 = sin(pkin(6));
	t180 = t169 * t165;
	t171 = cos(pkin(6));
	t172 = sin(qJ(2));
	t179 = t171 * t172;
	t173 = cos(qJ(2));
	t178 = t171 * t173;
	t177 = qJD(2) * t165;
	t176 = qJD(2) * t169;
	t168 = sin(pkin(11));
	t170 = cos(pkin(11));
	t175 = t168 * t173 + t170 * t179;
	t174 = -t168 * t179 + t170 * t173;
	t166 = cos(t167);
	t164 = t173 * t165 * t176 + (t166 * t169 * t172 + t165 * t171) * qJD(3);
	t163 = (t174 * t166 + t168 * t180) * qJD(3) + (-t168 * t178 - t170 * t172) * t177;
	t162 = (t175 * t166 - t170 * t180) * qJD(3) + (-t168 * t172 + t170 * t178) * t177;
	t1 = [0, 0, t174 * qJD(2), 0, t163, t163; 0, 0, t175 * qJD(2), 0, t162, t162; 0, 0, t172 * t176, 0, t164, t164;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end