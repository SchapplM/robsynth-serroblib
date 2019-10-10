% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRPRR1
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
% Datum: 2019-10-09 22:25
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRPRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:35
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:35
	% EndTime: 2019-10-09 22:25:36
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
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
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
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t116 = sin(qJ(2));
	t118 = cos(pkin(6)) * t116;
	t117 = cos(qJ(2));
	t114 = cos(pkin(11));
	t113 = sin(pkin(11));
	t112 = sin(pkin(6)) * qJD(2) * t116;
	t111 = (-t113 * t118 + t114 * t117) * qJD(2);
	t110 = (t113 * t117 + t114 * t118) * qJD(2);
	t1 = [0, 0, t111, 0, t111, 0; 0, 0, t110, 0, t110, 0; 0, 0, t112, 0, t112, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 22:25:36
	% EndTime: 2019-10-09 22:25:36
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (38->12), mult. (60->30), div. (0->0), fcn. (62->8), ass. (0->21)
	t174 = qJ(3) + pkin(12) + qJ(5);
	t172 = sin(t174);
	t177 = sin(pkin(6));
	t188 = t177 * t172;
	t179 = cos(pkin(6));
	t180 = sin(qJ(2));
	t187 = t179 * t180;
	t181 = cos(qJ(2));
	t186 = t179 * t181;
	t185 = qJD(2) * t172;
	t184 = qJD(2) * t177;
	t176 = sin(pkin(11));
	t178 = cos(pkin(11));
	t183 = t176 * t181 + t178 * t187;
	t182 = -t176 * t187 + t178 * t181;
	t175 = qJD(3) + qJD(5);
	t173 = cos(t174);
	t171 = t180 * t184;
	t170 = t182 * qJD(2);
	t169 = t183 * qJD(2);
	t1 = [0, 0, t170, 0, t170, (t173 * t182 + t176 * t188) * t175 + (-t176 * t186 - t178 * t180) * t185; 0, 0, t169, 0, t169, (t173 * t183 - t178 * t188) * t175 + (-t176 * t180 + t178 * t186) * t185; 0, 0, t171, 0, t171, t177 * t180 * t175 * t173 + (t175 * t179 + t181 * t184) * t172;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end