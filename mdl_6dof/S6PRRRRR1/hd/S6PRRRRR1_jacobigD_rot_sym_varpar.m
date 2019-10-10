% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6PRRRRR1
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
% Datum: 2019-10-09 23:13
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRR1_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6PRRRRR1_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobigD_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:44
	% EndTime: 2019-10-09 23:13:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
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
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (4->2), mult. (24->9), div. (0->0), fcn. (24->6), ass. (0->9)
	t111 = sin(qJ(2));
	t113 = cos(pkin(6)) * t111;
	t112 = cos(qJ(2));
	t109 = cos(pkin(12));
	t108 = sin(pkin(12));
	t107 = sin(pkin(6)) * qJD(2) * t111;
	t106 = (-t108 * t113 + t109 * t112) * qJD(2);
	t105 = (t108 * t112 + t109 * t113) * qJD(2);
	t1 = [0, 0, t106, t106, 0, 0; 0, 0, t105, t105, 0, 0; 0, 0, t107, t107, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->2), mult. (36->9), div. (0->0), fcn. (36->6), ass. (0->9)
	t115 = sin(qJ(2));
	t117 = cos(pkin(6)) * t115;
	t116 = cos(qJ(2));
	t113 = cos(pkin(12));
	t112 = sin(pkin(12));
	t111 = sin(pkin(6)) * qJD(2) * t115;
	t110 = (-t112 * t117 + t113 * t116) * qJD(2);
	t109 = (t112 * t116 + t113 * t117) * qJD(2);
	t1 = [0, 0, t110, t110, t110, 0; 0, 0, t109, t109, t109, 0; 0, 0, t111, t111, t111, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 23:13:45
	% EndTime: 2019-10-09 23:13:45
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (46->12), mult. (72->30), div. (0->0), fcn. (74->8), ass. (0->21)
	t178 = qJ(3) + qJ(4) + qJ(5);
	t175 = sin(t178);
	t180 = sin(pkin(6));
	t191 = t180 * t175;
	t182 = cos(pkin(6));
	t183 = sin(qJ(2));
	t190 = t182 * t183;
	t184 = cos(qJ(2));
	t189 = t182 * t184;
	t188 = qJD(2) * t175;
	t187 = qJD(2) * t180;
	t179 = sin(pkin(12));
	t181 = cos(pkin(12));
	t186 = t179 * t184 + t181 * t190;
	t185 = -t179 * t190 + t181 * t184;
	t177 = qJD(3) + qJD(4) + qJD(5);
	t176 = cos(t178);
	t174 = t183 * t187;
	t173 = t185 * qJD(2);
	t172 = t186 * qJD(2);
	t1 = [0, 0, t173, t173, t173, (t185 * t176 + t179 * t191) * t177 + (-t179 * t189 - t181 * t183) * t188; 0, 0, t172, t172, t172, (t186 * t176 - t181 * t191) * t177 + (-t179 * t183 + t181 * t189) * t188; 0, 0, t174, t174, t174, t180 * t183 * t177 * t176 + (t177 * t182 + t184 * t187) * t175;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end