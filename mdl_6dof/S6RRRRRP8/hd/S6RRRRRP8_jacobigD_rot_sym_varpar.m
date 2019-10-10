% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:07
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP8_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP8_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP8_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP8_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
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
	% StartTime: 2019-10-10 13:07:41
	% EndTime: 2019-10-10 13:07:41
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (12->6), mult. (48->16), div. (0->0), fcn. (48->6), ass. (0->15)
	t128 = sin(qJ(2));
	t129 = sin(qJ(1));
	t136 = t128 * t129;
	t131 = cos(qJ(1));
	t135 = t128 * t131;
	t130 = cos(qJ(2));
	t134 = t129 * t130;
	t133 = t130 * t131;
	t126 = sin(pkin(6));
	t132 = qJD(1) * t126;
	t127 = cos(pkin(6));
	t125 = t126 * qJD(2) * t128;
	t124 = (t127 * t135 + t134) * qJD(2) + (t127 * t134 + t135) * qJD(1);
	t123 = (-t127 * t136 + t133) * qJD(2) + (t127 * t133 - t136) * qJD(1);
	t1 = [0, t131 * t132, t123, t123, 0, 0; 0, t129 * t132, t124, t124, 0, 0; 0, 0, t125, t125, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:42
	% EndTime: 2019-10-10 13:07:42
	% DurationCPUTime: 0.16s
	% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
	t197 = sin(pkin(6));
	t200 = sin(qJ(1));
	t215 = t197 * t200;
	t202 = cos(qJ(1));
	t214 = t197 * t202;
	t199 = sin(qJ(2));
	t213 = t199 * t200;
	t212 = t199 * t202;
	t201 = cos(qJ(2));
	t211 = t200 * t201;
	t210 = t202 * t201;
	t209 = qJD(1) * t197;
	t196 = qJ(3) + qJ(4);
	t193 = sin(t196);
	t208 = qJD(2) * t193;
	t207 = qJD(2) * t197;
	t198 = cos(pkin(6));
	t206 = t198 * t210 - t213;
	t205 = t198 * t211 + t212;
	t204 = t198 * t212 + t211;
	t203 = -t198 * t213 + t210;
	t195 = qJD(3) + qJD(4);
	t194 = cos(t196);
	t192 = t199 * t207;
	t191 = qJD(1) * t205 + t204 * qJD(2);
	t190 = qJD(1) * t206 + t203 * qJD(2);
	t1 = [0, t202 * t209, t190, t190, (t193 * t215 + t194 * t203) * t195 - t205 * t208 + (-t193 * t204 - t194 * t214) * qJD(1), 0; 0, t200 * t209, t191, t191, (-t193 * t214 + t194 * t204) * t195 + t206 * t208 + (t193 * t203 - t194 * t215) * qJD(1), 0; 0, 0, t192, t192, t197 * t199 * t195 * t194 + (t195 * t198 + t201 * t207) * t193, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:07:42
	% EndTime: 2019-10-10 13:07:42
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (45->18), mult. (100->40), div. (0->0), fcn. (102->8), ass. (0->27)
	t240 = sin(pkin(6));
	t243 = sin(qJ(1));
	t258 = t240 * t243;
	t245 = cos(qJ(1));
	t257 = t240 * t245;
	t242 = sin(qJ(2));
	t256 = t242 * t243;
	t255 = t242 * t245;
	t244 = cos(qJ(2));
	t254 = t243 * t244;
	t253 = t245 * t244;
	t252 = qJD(1) * t240;
	t239 = qJ(3) + qJ(4);
	t236 = sin(t239);
	t251 = qJD(2) * t236;
	t250 = qJD(2) * t240;
	t241 = cos(pkin(6));
	t249 = t241 * t253 - t256;
	t248 = t241 * t254 + t255;
	t247 = t241 * t255 + t254;
	t246 = -t241 * t256 + t253;
	t238 = qJD(3) + qJD(4);
	t237 = cos(t239);
	t235 = t242 * t250;
	t234 = t248 * qJD(1) + t247 * qJD(2);
	t233 = t249 * qJD(1) + t246 * qJD(2);
	t1 = [0, t245 * t252, t233, t233, (t236 * t258 + t246 * t237) * t238 - t248 * t251 + (-t247 * t236 - t237 * t257) * qJD(1), 0; 0, t243 * t252, t234, t234, (-t236 * t257 + t247 * t237) * t238 + t249 * t251 + (t246 * t236 - t237 * t258) * qJD(1), 0; 0, 0, t235, t235, t240 * t242 * t238 * t237 + (t238 * t241 + t244 * t250) * t236, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end