% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 10:35
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRRP5_jacobigD_rot_sym_varpar(qJ, qJD, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_jacobigD_rot_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRPRRP5_jacobigD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_jacobigD_rot_sym_varpar: pkin has to be [11x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobigD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobigD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobigD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t59 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t59, 0, 0, 0, 0; 0, sin(qJ(1)) * t59, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobigD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.05s
	% Computational Cost: add. (0->0), mult. (4->3), div. (0->0), fcn. (4->3), ass. (0->2)
	t79 = qJD(1) * sin(pkin(6));
	t1 = [0, cos(qJ(1)) * t79, 0, 0, 0, 0; 0, sin(qJ(1)) * t79, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	JgD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobigD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:23
	% EndTime: 2019-10-10 10:35:23
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (15->8), mult. (54->22), div. (0->0), fcn. (58->8), ass. (0->16)
	t133 = sin(pkin(11));
	t135 = cos(pkin(11));
	t137 = sin(qJ(2));
	t139 = cos(qJ(2));
	t142 = t139 * t133 + t137 * t135;
	t144 = qJD(2) * t142;
	t134 = sin(pkin(6));
	t143 = qJD(1) * t134;
	t141 = t133 * t137 - t135 * t139;
	t140 = cos(qJ(1));
	t138 = sin(qJ(1));
	t136 = cos(pkin(6));
	t131 = t141 * qJD(2);
	t130 = t141 * t136;
	t129 = t136 * t144;
	t1 = [0, t140 * t143, 0, -t138 * t129 - t140 * t131 + (-t130 * t140 - t138 * t142) * qJD(1), 0, 0; 0, t138 * t143, 0, t140 * t129 - t138 * t131 + (-t130 * t138 + t140 * t142) * qJD(1), 0, 0; 0, 0, 0, t134 * t144, 0, 0;];
	JgD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobigD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:24
	% EndTime: 2019-10-10 10:35:24
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
	t202 = sin(pkin(11));
	t204 = cos(pkin(11));
	t207 = sin(qJ(2));
	t210 = cos(qJ(2));
	t213 = t207 * t202 - t210 * t204;
	t220 = t213 * qJD(2);
	t214 = t210 * t202 + t207 * t204;
	t199 = t214 * qJD(2);
	t203 = sin(pkin(6));
	t208 = sin(qJ(1));
	t219 = t203 * t208;
	t211 = cos(qJ(1));
	t218 = t203 * t211;
	t217 = qJD(1) * t203;
	t205 = cos(pkin(6));
	t197 = t214 * t205;
	t216 = t211 * t197 - t208 * t213;
	t215 = -t208 * t197 - t211 * t213;
	t209 = cos(qJ(4));
	t206 = sin(qJ(4));
	t196 = t213 * t205;
	t195 = t205 * t220;
	t194 = t205 * t199;
	t1 = [0, t211 * t217, 0, -t208 * t194 - t211 * t220 + (-t196 * t211 - t208 * t214) * qJD(1), (t208 * t195 - t211 * t199) * t206 + (t206 * t219 + t215 * t209) * qJD(4) + (-t216 * t206 - t209 * t218) * qJD(1), 0; 0, t208 * t217, 0, t211 * t194 - t208 * t220 + (-t196 * t208 + t211 * t214) * qJD(1), (-t211 * t195 - t208 * t199) * t206 + (-t206 * t218 + t216 * t209) * qJD(4) + (t215 * t206 - t209 * t219) * qJD(1), 0; 0, 0, 0, t203 * t199, t205 * qJD(4) * t206 + (t214 * qJD(4) * t209 - t206 * t220) * t203, 0;];
	JgD_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobigD_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 10:35:24
	% EndTime: 2019-10-10 10:35:24
	% DurationCPUTime: 0.25s
	% Computational Cost: add. (45->22), mult. (147->54), div. (0->0), fcn. (162->10), ass. (0->24)
	t204 = sin(pkin(11));
	t206 = cos(pkin(11));
	t209 = sin(qJ(2));
	t212 = cos(qJ(2));
	t215 = t209 * t204 - t212 * t206;
	t222 = t215 * qJD(2);
	t216 = t212 * t204 + t209 * t206;
	t201 = t216 * qJD(2);
	t205 = sin(pkin(6));
	t210 = sin(qJ(1));
	t221 = t205 * t210;
	t213 = cos(qJ(1));
	t220 = t205 * t213;
	t219 = qJD(1) * t205;
	t207 = cos(pkin(6));
	t199 = t216 * t207;
	t218 = t213 * t199 - t210 * t215;
	t217 = -t210 * t199 - t213 * t215;
	t211 = cos(qJ(4));
	t208 = sin(qJ(4));
	t198 = t215 * t207;
	t197 = t207 * t222;
	t196 = t207 * t201;
	t1 = [0, t213 * t219, 0, -t210 * t196 - t213 * t222 + (-t198 * t213 - t210 * t216) * qJD(1), (t210 * t197 - t213 * t201) * t208 + (t208 * t221 + t217 * t211) * qJD(4) + (-t218 * t208 - t211 * t220) * qJD(1), 0; 0, t210 * t219, 0, t213 * t196 - t210 * t222 + (-t198 * t210 + t213 * t216) * qJD(1), (-t213 * t197 - t210 * t201) * t208 + (-t208 * t220 + t218 * t211) * qJD(4) + (t217 * t208 - t211 * t221) * qJD(1), 0; 0, 0, 0, t205 * t201, t207 * qJD(4) * t208 + (t216 * qJD(4) * t211 - t208 * t222) * t205, 0;];
	JgD_rot = t1;
else
	JgD_rot=NaN(3,6);
end