% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S5PRRPR7
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
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% JgD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S5PRRPR7_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_jacobigD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_jacobigD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRPR7_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:31
	% EndTime: 2019-12-05 16:38:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t82 = sin(qJ(2));
	t84 = cos(pkin(5)) * t82;
	t83 = cos(qJ(2));
	t80 = cos(pkin(9));
	t79 = sin(pkin(9));
	t1 = [0, 0, (-t79 * t84 + t80 * t83) * qJD(2), 0, 0; 0, 0, (t79 * t83 + t80 * t84) * qJD(2), 0, 0; 0, 0, sin(pkin(5)) * qJD(2) * t82, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (2->2), mult. (12->9), div. (0->0), fcn. (12->6), ass. (0->6)
	t116 = sin(qJ(2));
	t118 = cos(pkin(5)) * t116;
	t117 = cos(qJ(2));
	t114 = cos(pkin(9));
	t113 = sin(pkin(9));
	t1 = [0, 0, (-t113 * t118 + t114 * t117) * qJD(2), 0, 0; 0, 0, (t113 * t117 + t114 * t118) * qJD(2), 0, 0; 0, 0, sin(pkin(5)) * qJD(2) * t116, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 16:38:32
	% EndTime: 2019-12-05 16:38:32
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (17->13), mult. (66->37), div. (0->0), fcn. (68->10), ass. (0->20)
	t176 = sin(pkin(5));
	t182 = cos(qJ(3));
	t190 = t176 * t182;
	t179 = cos(pkin(5));
	t181 = sin(qJ(2));
	t189 = t179 * t181;
	t183 = cos(qJ(2));
	t188 = t179 * t183;
	t187 = qJD(2) * t181;
	t186 = qJD(2) * t182;
	t175 = sin(pkin(9));
	t178 = cos(pkin(9));
	t185 = t175 * t183 + t178 * t189;
	t184 = -t175 * t189 + t178 * t183;
	t180 = sin(qJ(3));
	t177 = cos(pkin(10));
	t174 = sin(pkin(10));
	t173 = t184 * qJD(2);
	t172 = t185 * qJD(2);
	t1 = [0, 0, t173, 0, ((t175 * t190 - t184 * t180) * qJD(3) + (-t175 * t188 - t178 * t181) * t186) * t174 - t173 * t177; 0, 0, t172, 0, ((-t178 * t190 - t185 * t180) * qJD(3) + (-t175 * t181 + t178 * t188) * t186) * t174 - t172 * t177; 0, 0, t176 * t187, 0, t179 * qJD(3) * t182 * t174 + ((-qJD(3) * t180 * t181 + t183 * t186) * t174 - t177 * t187) * t176;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,5);
end