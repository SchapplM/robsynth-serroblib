% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRRP12
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: bsp_3T1R_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 13:16
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRRP12_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRRP12_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:47
	% EndTime: 2019-10-10 13:15:47
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t94 = cos(pkin(6));
	t97 = cos(qJ(2));
	t101 = t94 * t97;
	t92 = sin(pkin(6));
	t96 = sin(qJ(1));
	t100 = t96 * t92;
	t98 = cos(qJ(1));
	t99 = t98 * t92;
	t95 = sin(qJ(2));
	t93 = cos(pkin(7));
	t91 = sin(pkin(7));
	t1 = [0, t100, -(-t96 * t101 - t98 * t95) * t91 + t93 * t100, 0, 0, 0; 0, -t99, -(t98 * t101 - t96 * t95) * t91 - t93 * t99, 0, 0, 0; 1, t94, -t92 * t97 * t91 + t94 * t93, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:48
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
	t140 = sin(pkin(6));
	t145 = sin(qJ(1));
	t154 = t145 * t140;
	t144 = sin(qJ(2));
	t153 = t145 * t144;
	t147 = cos(qJ(2));
	t152 = t145 * t147;
	t148 = cos(qJ(1));
	t151 = t148 * t140;
	t150 = t148 * t144;
	t149 = t148 * t147;
	t146 = cos(qJ(3));
	t143 = sin(qJ(3));
	t142 = cos(pkin(6));
	t141 = cos(pkin(7));
	t139 = sin(pkin(7));
	t138 = -t142 * t152 - t150;
	t137 = t142 * t149 - t153;
	t1 = [0, t154, -t138 * t139 + t141 * t154, (-t142 * t153 + t149) * t143 + (-t138 * t141 - t139 * t154) * t146, 0, 0; 0, -t151, -t137 * t139 - t141 * t151, (t142 * t150 + t152) * t143 + (-t137 * t141 + t139 * t151) * t146, 0, 0; 1, t142, -t140 * t147 * t139 + t142 * t141, -t142 * t139 * t146 + (-t141 * t146 * t147 + t143 * t144) * t140, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:48
	% EndTime: 2019-10-10 13:15:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
	t187 = sin(pkin(7));
	t190 = cos(pkin(6));
	t208 = t187 * t190;
	t189 = cos(pkin(7));
	t197 = cos(qJ(2));
	t207 = t189 * t197;
	t188 = sin(pkin(6));
	t194 = sin(qJ(1));
	t206 = t194 * t188;
	t193 = sin(qJ(2));
	t205 = t194 * t193;
	t204 = t194 * t197;
	t198 = cos(qJ(1));
	t203 = t198 * t188;
	t202 = t198 * t193;
	t201 = t198 * t197;
	t183 = t190 * t201 - t205;
	t200 = -t183 * t189 + t187 * t203;
	t185 = -t190 * t204 - t202;
	t199 = t185 * t189 + t187 * t206;
	t196 = cos(qJ(3));
	t195 = cos(qJ(4));
	t192 = sin(qJ(3));
	t191 = sin(qJ(4));
	t186 = -t190 * t205 + t201;
	t184 = t190 * t202 + t204;
	t182 = -t188 * t197 * t187 + t190 * t189;
	t181 = -t185 * t187 + t189 * t206;
	t180 = -t183 * t187 - t189 * t203;
	t1 = [0, t206, t181, t186 * t192 - t199 * t196, (t186 * t196 + t199 * t192) * t191 - t181 * t195, 0; 0, -t203, t180, t184 * t192 + t200 * t196, (t184 * t196 - t200 * t192) * t191 - t180 * t195, 0; 1, t190, t182, -t196 * t208 + (t192 * t193 - t196 * t207) * t188, (t192 * t208 + (t192 * t207 + t193 * t196) * t188) * t191 - t182 * t195, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 13:15:49
	% EndTime: 2019-10-10 13:15:49
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
	t226 = sin(pkin(7));
	t229 = cos(pkin(6));
	t247 = t226 * t229;
	t228 = cos(pkin(7));
	t236 = cos(qJ(2));
	t246 = t228 * t236;
	t227 = sin(pkin(6));
	t233 = sin(qJ(1));
	t245 = t233 * t227;
	t232 = sin(qJ(2));
	t244 = t233 * t232;
	t243 = t233 * t236;
	t237 = cos(qJ(1));
	t242 = t237 * t227;
	t241 = t237 * t232;
	t240 = t237 * t236;
	t222 = t229 * t240 - t244;
	t239 = -t222 * t228 + t226 * t242;
	t224 = -t229 * t243 - t241;
	t238 = t224 * t228 + t226 * t245;
	t235 = cos(qJ(3));
	t234 = cos(qJ(4));
	t231 = sin(qJ(3));
	t230 = sin(qJ(4));
	t225 = -t229 * t244 + t240;
	t223 = t229 * t241 + t243;
	t221 = -t227 * t236 * t226 + t229 * t228;
	t220 = -t224 * t226 + t228 * t245;
	t219 = -t222 * t226 - t228 * t242;
	t1 = [0, t245, t220, t225 * t231 - t238 * t235, (t225 * t235 + t238 * t231) * t230 - t220 * t234, 0; 0, -t242, t219, t223 * t231 + t239 * t235, (t223 * t235 - t239 * t231) * t230 - t219 * t234, 0; 1, t229, t221, -t235 * t247 + (t231 * t232 - t235 * t246) * t227, (t231 * t247 + (t231 * t246 + t232 * t235) * t227) * t230 - t221 * t234, 0;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end