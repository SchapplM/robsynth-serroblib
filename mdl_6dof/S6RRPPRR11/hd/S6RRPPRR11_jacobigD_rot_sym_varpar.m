% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta4]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:53
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR11_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR11_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR11_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR11_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR11_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:50
	% EndTime: 2019-10-10 09:53:50
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t71 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t71, 0, 0, 0, 0; 0, sin(qJ(1)) * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t82 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t82, 0, 0, 0, 0; 0, sin(qJ(1)) * t82, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t110 = sin(qJ(2));
	t111 = sin(qJ(1));
	t118 = t110 * t111;
	t113 = cos(qJ(1));
	t117 = t110 * t113;
	t112 = cos(qJ(2));
	t116 = t111 * t112;
	t115 = t112 * t113;
	t108 = sin(pkin(6));
	t114 = qJD(1) * t108;
	t109 = cos(pkin(6));
	t1 = [0, t113 * t114, 0, 0, (-t109 * t116 - t117) * qJD(2) + (-t109 * t117 - t116) * qJD(1), 0; 0, t111 * t114, 0, 0, (t109 * t115 - t118) * qJD(2) + (-t109 * t118 + t115) * qJD(1), 0; 0, 0, 0, 0, t108 * qJD(2) * t112, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:53:51
	% EndTime: 2019-10-10 09:53:51
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
	t166 = sin(pkin(6));
	t169 = sin(qJ(1));
	t184 = t166 * t169;
	t171 = cos(qJ(1));
	t183 = t166 * t171;
	t168 = sin(qJ(2));
	t182 = t169 * t168;
	t170 = cos(qJ(2));
	t181 = t169 * t170;
	t180 = t170 * t171;
	t179 = t171 * t168;
	t178 = qJD(1) * t166;
	t165 = pkin(11) + qJ(5);
	t164 = cos(t165);
	t177 = qJD(2) * t164;
	t176 = qJD(2) * t166;
	t167 = cos(pkin(6));
	t175 = t167 * t180 - t182;
	t174 = t167 * t181 + t179;
	t173 = t167 * t179 + t181;
	t172 = -t167 * t182 + t180;
	t163 = sin(t165);
	t1 = [0, t171 * t178, 0, 0, -t173 * qJD(1) - t174 * qJD(2), (t174 * t163 + t164 * t184) * qJD(5) - t172 * t177 + (t163 * t183 - t175 * t164) * qJD(1); 0, t169 * t178, 0, 0, t172 * qJD(1) + t175 * qJD(2), (-t175 * t163 - t164 * t183) * qJD(5) - t173 * t177 + (t163 * t184 - t174 * t164) * qJD(1); 0, 0, 0, 0, t170 * t176, -t168 * t164 * t176 + (-t163 * t166 * t170 + t164 * t167) * qJD(5);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end