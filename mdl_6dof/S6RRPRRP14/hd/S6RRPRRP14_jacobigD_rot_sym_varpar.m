% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:50
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP14_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP14_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP14_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP14_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP14_jacobigD_rot_sym_varpar: pkin has to be [10x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t71 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t71, 0, 0, 0, 0; 0, sin(qJ(1)) * t71, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
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
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:26
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t155 = sin(pkin(6));
	t159 = sin(qJ(1));
	t175 = t155 * t159;
	t161 = cos(qJ(2));
	t174 = t155 * t161;
	t162 = cos(qJ(1));
	t173 = t155 * t162;
	t158 = sin(qJ(2));
	t172 = t159 * t158;
	t171 = t159 * t161;
	t170 = t161 * t162;
	t169 = t162 * t158;
	t168 = qJD(1) * t155;
	t160 = cos(qJ(4));
	t167 = qJD(2) * t160;
	t156 = cos(pkin(6));
	t166 = t156 * t170 - t172;
	t165 = t156 * t171 + t169;
	t164 = t156 * t169 + t171;
	t163 = -t156 * t172 + t170;
	t157 = sin(qJ(4));
	t1 = [0, t162 * t168, 0, -t164 * qJD(1) - t165 * qJD(2), (t165 * t157 + t160 * t175) * qJD(4) - t163 * t167 + (t157 * t173 - t166 * t160) * qJD(1), 0; 0, t159 * t168, 0, t163 * qJD(1) + t166 * qJD(2), (-t166 * t157 - t160 * t173) * qJD(4) - t164 * t167 + (t157 * t175 - t165 * t160) * qJD(1), 0; 0, 0, 0, qJD(2) * t174, -t155 * t158 * t167 + (t156 * t160 - t157 * t174) * qJD(4), 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:50:26
	% EndTime: 2019-10-10 10:50:27
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
	t192 = sin(pkin(6));
	t196 = sin(qJ(1));
	t212 = t192 * t196;
	t198 = cos(qJ(2));
	t211 = t192 * t198;
	t199 = cos(qJ(1));
	t210 = t192 * t199;
	t195 = sin(qJ(2));
	t209 = t196 * t195;
	t208 = t196 * t198;
	t207 = t198 * t199;
	t206 = t199 * t195;
	t205 = qJD(1) * t192;
	t197 = cos(qJ(4));
	t204 = qJD(2) * t197;
	t193 = cos(pkin(6));
	t203 = t193 * t207 - t209;
	t202 = t193 * t208 + t206;
	t201 = t193 * t206 + t208;
	t200 = -t193 * t209 + t207;
	t194 = sin(qJ(4));
	t1 = [0, t199 * t205, 0, -t201 * qJD(1) - t202 * qJD(2), (t202 * t194 + t197 * t212) * qJD(4) - t200 * t204 + (t194 * t210 - t203 * t197) * qJD(1), 0; 0, t196 * t205, 0, t200 * qJD(1) + t203 * qJD(2), (-t203 * t194 - t197 * t210) * qJD(4) - t201 * t204 + (t194 * t212 - t202 * t197) * qJD(1), 0; 0, 0, 0, qJD(2) * t211, -t192 * t195 * t204 + (t193 * t197 - t194 * t211) * qJD(4), 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end