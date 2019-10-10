% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPPRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 09:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPPRR9_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR9_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR9_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPPRR9_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR9_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t71 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t71, 0, 0, 0, 0; 0, sin(qJ(1)) * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t69 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t69, 0, 0, 0, 0; 0, sin(qJ(1)) * t69, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (7->7), mult. (26->17), div. (0->0), fcn. (26->6), ass. (0->12)
	t106 = sin(qJ(2));
	t107 = sin(qJ(1));
	t114 = t106 * t107;
	t109 = cos(qJ(1));
	t113 = t106 * t109;
	t108 = cos(qJ(2));
	t112 = t107 * t108;
	t111 = t108 * t109;
	t104 = sin(pkin(6));
	t110 = qJD(1) * t104;
	t105 = cos(pkin(6));
	t1 = [0, t109 * t110, 0, 0, (t105 * t114 - t111) * qJD(2) + (-t105 * t111 + t114) * qJD(1), 0; 0, t107 * t110, 0, 0, (-t105 * t113 - t112) * qJD(2) + (-t105 * t112 - t113) * qJD(1), 0; 0, 0, 0, 0, -t104 * qJD(2) * t106, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 09:50:17
	% EndTime: 2019-10-10 09:50:17
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (23->17), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t157 = sin(pkin(6));
	t160 = sin(qJ(2));
	t177 = t157 * t160;
	t161 = sin(qJ(1));
	t176 = t157 * t161;
	t164 = cos(qJ(1));
	t175 = t157 * t164;
	t174 = t160 * t161;
	t173 = t160 * t164;
	t163 = cos(qJ(2));
	t172 = t161 * t163;
	t171 = t164 * t163;
	t170 = qJD(1) * t157;
	t162 = cos(qJ(5));
	t169 = qJD(2) * t162;
	t158 = cos(pkin(6));
	t168 = t158 * t171 - t174;
	t167 = -t158 * t172 - t173;
	t166 = t158 * t173 + t172;
	t165 = t158 * t174 - t171;
	t159 = sin(qJ(5));
	t1 = [0, t164 * t170, 0, 0, -t168 * qJD(1) + t165 * qJD(2), (-t165 * t159 + t162 * t176) * qJD(5) - t167 * t169 + (t159 * t175 + t166 * t162) * qJD(1); 0, t161 * t170, 0, 0, t167 * qJD(1) - t166 * qJD(2), (t166 * t159 - t162 * t175) * qJD(5) - t168 * t169 + (t159 * t176 + t165 * t162) * qJD(1); 0, 0, 0, 0, -qJD(2) * t177, -t157 * t163 * t169 + (t158 * t162 + t159 * t177) * qJD(5);];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end