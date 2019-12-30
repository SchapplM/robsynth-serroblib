% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 20:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5RRRPR10_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR10_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRPR10_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:35
	% EndTime: 2019-12-29 20:13:35
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:37
	% EndTime: 2019-12-29 20:13:37
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:37
	% EndTime: 2019-12-29 20:13:37
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(5));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:37
	% EndTime: 2019-12-29 20:13:37
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t103 = sin(qJ(2));
	t104 = sin(qJ(1));
	t111 = t103 * t104;
	t106 = cos(qJ(1));
	t110 = t103 * t106;
	t105 = cos(qJ(2));
	t109 = t104 * t105;
	t108 = t105 * t106;
	t101 = sin(pkin(5));
	t107 = qJD(1) * t101;
	t102 = cos(pkin(5));
	t1 = [0, t106 * t107, (-t102 * t111 + t108) * qJD(2) + (t102 * t108 - t111) * qJD(1), 0, 0; 0, t104 * t107, (t102 * t110 + t109) * qJD(2) + (t102 * t109 + t110) * qJD(1), 0, 0; 0, 0, t101 * qJD(2) * t103, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:37
	% EndTime: 2019-12-29 20:13:37
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t109 = sin(qJ(2));
	t110 = sin(qJ(1));
	t117 = t109 * t110;
	t112 = cos(qJ(1));
	t116 = t109 * t112;
	t111 = cos(qJ(2));
	t115 = t110 * t111;
	t114 = t111 * t112;
	t107 = sin(pkin(5));
	t113 = qJD(1) * t107;
	t108 = cos(pkin(5));
	t1 = [0, t112 * t113, (-t108 * t117 + t114) * qJD(2) + (t108 * t114 - t117) * qJD(1), 0, 0; 0, t110 * t113, (t108 * t116 + t115) * qJD(2) + (t108 * t115 + t116) * qJD(1), 0, 0; 0, 0, t107 * qJD(2) * t109, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 20:13:36
	% EndTime: 2019-12-29 20:13:36
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t168 = sin(pkin(5));
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
	t167 = qJ(3) + pkin(10);
	t165 = sin(t167);
	t179 = qJD(2) * t165;
	t178 = qJD(2) * t168;
	t169 = cos(pkin(5));
	t177 = t169 * t181 - t184;
	t176 = t169 * t182 + t183;
	t175 = t169 * t183 + t182;
	t174 = -t169 * t184 + t181;
	t166 = cos(t167);
	t1 = [0, t173 * t180, t177 * qJD(1) + t174 * qJD(2), 0, (t165 * t186 + t174 * t166) * qJD(3) - t176 * t179 + (-t175 * t165 - t166 * t185) * qJD(1); 0, t171 * t180, t176 * qJD(1) + t175 * qJD(2), 0, (-t165 * t185 + t175 * t166) * qJD(3) + t177 * t179 + (t174 * t165 - t166 * t186) * qJD(1); 0, 0, t170 * t178, 0, t172 * t165 * t178 + (t166 * t168 * t170 + t165 * t169) * qJD(3);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,5);
end