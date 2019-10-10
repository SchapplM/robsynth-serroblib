% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR15
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:18
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR15_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR15_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobig_rot_sym_varpar: pkin has to be [12x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t123 = cos(pkin(6));
	t126 = cos(qJ(2));
	t130 = t123 * t126;
	t121 = sin(pkin(6));
	t125 = sin(qJ(1));
	t129 = t125 * t121;
	t127 = cos(qJ(1));
	t128 = t127 * t121;
	t124 = sin(qJ(2));
	t122 = cos(pkin(7));
	t120 = sin(pkin(7));
	t1 = [0, t129, -(-t127 * t124 - t125 * t130) * t120 + t122 * t129, 0, 0, 0; 0, -t128, -(-t125 * t124 + t127 * t130) * t120 - t122 * t128, 0, 0, 0; 1, t123, -t121 * t126 * t120 + t123 * t122, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:23
	% EndTime: 2019-10-10 12:18:23
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
	t1 = [0, t154, -t138 * t139 + t141 * t154, 0, (-t142 * t153 + t149) * t146 + (t138 * t141 + t139 * t154) * t143, 0; 0, -t151, -t137 * t139 - t141 * t151, 0, (t142 * t150 + t152) * t146 + (t137 * t141 - t139 * t151) * t143, 0; 1, t142, -t140 * t147 * t139 + t142 * t141, 0, t142 * t139 * t143 + (t141 * t143 * t147 + t144 * t146) * t140, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:18:24
	% EndTime: 2019-10-10 12:18:24
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
	t188 = sin(pkin(7));
	t191 = cos(pkin(6));
	t209 = t188 * t191;
	t190 = cos(pkin(7));
	t198 = cos(qJ(2));
	t208 = t190 * t198;
	t189 = sin(pkin(6));
	t195 = sin(qJ(1));
	t207 = t195 * t189;
	t194 = sin(qJ(2));
	t206 = t195 * t194;
	t205 = t195 * t198;
	t199 = cos(qJ(1));
	t204 = t199 * t189;
	t203 = t199 * t194;
	t202 = t199 * t198;
	t184 = t191 * t202 - t206;
	t201 = -t184 * t190 + t188 * t204;
	t186 = -t191 * t205 - t203;
	t200 = t186 * t190 + t188 * t207;
	t197 = cos(qJ(3));
	t196 = cos(qJ(5));
	t193 = sin(qJ(3));
	t192 = sin(qJ(5));
	t187 = -t191 * t206 + t202;
	t185 = t191 * t203 + t205;
	t183 = -t189 * t198 * t188 + t191 * t190;
	t182 = -t186 * t188 + t190 * t207;
	t181 = -t184 * t188 - t190 * t204;
	t1 = [0, t207, t182, 0, t187 * t197 + t200 * t193, t182 * t192 - (t187 * t193 - t200 * t197) * t196; 0, -t204, t181, 0, t185 * t197 - t201 * t193, t181 * t192 - (t185 * t193 + t201 * t197) * t196; 1, t191, t183, 0, t193 * t209 + (t193 * t208 + t194 * t197) * t189, t183 * t192 - (-t197 * t209 + (t193 * t194 - t197 * t208) * t189) * t196;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end