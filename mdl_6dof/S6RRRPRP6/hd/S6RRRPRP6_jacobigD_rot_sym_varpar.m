% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 11:44
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRPRP6_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRP6_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
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
	% StartTime: 2019-10-10 11:43:59
	% EndTime: 2019-10-10 11:43:59
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t109 = sin(qJ(2));
	t110 = sin(qJ(1));
	t117 = t109 * t110;
	t112 = cos(qJ(1));
	t116 = t109 * t112;
	t111 = cos(qJ(2));
	t115 = t110 * t111;
	t114 = t111 * t112;
	t107 = sin(pkin(6));
	t113 = qJD(1) * t107;
	t108 = cos(pkin(6));
	t1 = [0, t112 * t113, (-t108 * t117 + t114) * qJD(2) + (t108 * t114 - t117) * qJD(1), 0, 0, 0; 0, t110 * t113, (t108 * t116 + t115) * qJD(2) + (t108 * t115 + t116) * qJD(1), 0, 0, 0; 0, 0, t107 * qJD(2) * t109, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:44:00
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t168 = sin(pkin(6));
	t171 = sin(qJ(1));
	t186 = t168 * t171;
	t173 = cos(qJ(1));
	t185 = t168 * t173;
	t170 = sin(qJ(2));
	t184 = t170 * t171;
	t183 = t170 * t173;
	t172 = cos(qJ(2));
	t182 = t171 * t172;
	t181 = t173 * t172;
	t180 = qJD(1) * t168;
	t167 = qJ(3) + pkin(11);
	t165 = sin(t167);
	t179 = qJD(2) * t165;
	t178 = qJD(2) * t168;
	t169 = cos(pkin(6));
	t177 = t169 * t181 - t184;
	t176 = t169 * t182 + t183;
	t175 = t169 * t183 + t182;
	t174 = -t169 * t184 + t181;
	t166 = cos(t167);
	t1 = [0, t173 * t180, t177 * qJD(1) + t174 * qJD(2), 0, (t165 * t186 + t174 * t166) * qJD(3) - t176 * t179 + (-t175 * t165 - t166 * t185) * qJD(1), 0; 0, t171 * t180, t176 * qJD(1) + t175 * qJD(2), 0, (-t165 * t185 + t175 * t166) * qJD(3) + t177 * t179 + (t174 * t165 - t166 * t186) * qJD(1), 0; 0, 0, t170 * t178, 0, t172 * t165 * t178 + (t166 * t168 * t170 + t165 * t169) * qJD(3), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 11:44:00
	% EndTime: 2019-10-10 11:44:00
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t170 = sin(pkin(6));
	t173 = sin(qJ(1));
	t188 = t170 * t173;
	t175 = cos(qJ(1));
	t187 = t170 * t175;
	t172 = sin(qJ(2));
	t186 = t172 * t173;
	t185 = t172 * t175;
	t174 = cos(qJ(2));
	t184 = t173 * t174;
	t183 = t175 * t174;
	t182 = qJD(1) * t170;
	t169 = qJ(3) + pkin(11);
	t167 = sin(t169);
	t181 = qJD(2) * t167;
	t180 = qJD(2) * t170;
	t171 = cos(pkin(6));
	t179 = t171 * t183 - t186;
	t178 = t171 * t184 + t185;
	t177 = t171 * t185 + t184;
	t176 = -t171 * t186 + t183;
	t168 = cos(t169);
	t1 = [0, t175 * t182, t179 * qJD(1) + t176 * qJD(2), 0, (t167 * t188 + t176 * t168) * qJD(3) - t178 * t181 + (-t177 * t167 - t168 * t187) * qJD(1), 0; 0, t173 * t182, t178 * qJD(1) + t177 * qJD(2), 0, (-t167 * t187 + t177 * t168) * qJD(3) + t179 * t181 + (t176 * t167 - t168 * t188) * qJD(1), 0; 0, 0, t172 * t180, 0, t174 * t167 * t180 + (t168 * t170 * t172 + t167 * t171) * qJD(3), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end