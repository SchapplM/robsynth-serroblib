% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:24
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR12_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRPR12_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t71 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t71, 0, 0, 0, 0; 0, sin(qJ(1)) * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t105 = sin(qJ(2));
	t106 = sin(qJ(1));
	t113 = t105 * t106;
	t108 = cos(qJ(1));
	t112 = t105 * t108;
	t107 = cos(qJ(2));
	t111 = t106 * t107;
	t110 = t107 * t108;
	t103 = sin(pkin(6));
	t109 = qJD(1) * t103;
	t104 = cos(pkin(6));
	t1 = [0, t108 * t109, 0, (-t104 * t111 - t112) * qJD(2) + (-t104 * t112 - t111) * qJD(1), 0, 0; 0, t106 * t109, 0, (t104 * t110 - t113) * qJD(2) + (-t104 * t113 + t110) * qJD(1), 0, 0; 0, 0, 0, t103 * qJD(2) * t107, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t113 = sin(qJ(2));
	t114 = sin(qJ(1));
	t121 = t113 * t114;
	t116 = cos(qJ(1));
	t120 = t113 * t116;
	t115 = cos(qJ(2));
	t119 = t114 * t115;
	t118 = t115 * t116;
	t111 = sin(pkin(6));
	t117 = qJD(1) * t111;
	t112 = cos(pkin(6));
	t1 = [0, t116 * t117, 0, (-t112 * t119 - t120) * qJD(2) + (-t112 * t120 - t119) * qJD(1), 0, 0; 0, t114 * t117, 0, (t112 * t118 - t121) * qJD(2) + (-t112 * t121 + t118) * qJD(1), 0, 0; 0, 0, 0, t111 * qJD(2) * t115, 0, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:24:31
	% EndTime: 2019-10-10 10:24:31
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t167 = sin(pkin(6));
	t170 = sin(qJ(1));
	t185 = t167 * t170;
	t172 = cos(qJ(1));
	t184 = t167 * t172;
	t169 = sin(qJ(2));
	t183 = t170 * t169;
	t171 = cos(qJ(2));
	t182 = t170 * t171;
	t181 = t171 * t172;
	t180 = t172 * t169;
	t179 = qJD(1) * t167;
	t166 = qJ(4) + pkin(11);
	t165 = cos(t166);
	t178 = qJD(2) * t165;
	t177 = qJD(2) * t167;
	t168 = cos(pkin(6));
	t176 = t168 * t181 - t183;
	t175 = t168 * t182 + t180;
	t174 = t168 * t180 + t182;
	t173 = -t168 * t183 + t181;
	t164 = sin(t166);
	t1 = [0, t172 * t179, 0, -t174 * qJD(1) - t175 * qJD(2), 0, (t175 * t164 + t165 * t185) * qJD(4) - t173 * t178 + (t164 * t184 - t176 * t165) * qJD(1); 0, t170 * t179, 0, t173 * qJD(1) + t176 * qJD(2), 0, (-t176 * t164 - t165 * t184) * qJD(4) - t174 * t178 + (t164 * t185 - t175 * t165) * qJD(1); 0, 0, 0, t171 * t177, 0, -t169 * t165 * t177 + (-t164 * t167 * t171 + t165 * t168) * qJD(4);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end