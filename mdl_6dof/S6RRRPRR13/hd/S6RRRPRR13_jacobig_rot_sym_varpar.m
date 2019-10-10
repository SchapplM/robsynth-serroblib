% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR13
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:14
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR13_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR13_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:17
	% EndTime: 2019-10-10 12:14:17
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
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
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t136 = cos(pkin(6));
	t139 = cos(qJ(2));
	t143 = t136 * t139;
	t134 = sin(pkin(6));
	t138 = sin(qJ(1));
	t142 = t138 * t134;
	t140 = cos(qJ(1));
	t141 = t140 * t134;
	t137 = sin(qJ(2));
	t135 = cos(pkin(7));
	t133 = sin(pkin(7));
	t1 = [0, t142, -(-t140 * t137 - t138 * t143) * t133 + t135 * t142, 0, 0, 0; 0, -t141, -(-t138 * t137 + t140 * t143) * t133 - t135 * t141, 0, 0, 0; 1, t136, -t134 * t139 * t133 + t136 * t135, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:18
	% EndTime: 2019-10-10 12:14:18
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
	t147 = sin(pkin(6));
	t152 = sin(qJ(1));
	t161 = t152 * t147;
	t151 = sin(qJ(2));
	t160 = t152 * t151;
	t154 = cos(qJ(2));
	t159 = t152 * t154;
	t155 = cos(qJ(1));
	t158 = t155 * t147;
	t157 = t155 * t151;
	t156 = t155 * t154;
	t153 = cos(qJ(3));
	t150 = sin(qJ(3));
	t149 = cos(pkin(6));
	t148 = cos(pkin(7));
	t146 = sin(pkin(7));
	t145 = -t149 * t159 - t157;
	t144 = t149 * t156 - t160;
	t1 = [0, t161, -t145 * t146 + t148 * t161, 0, (-t149 * t160 + t156) * t150 + (-t145 * t148 - t146 * t161) * t153, 0; 0, -t158, -t144 * t146 - t148 * t158, 0, (t149 * t157 + t159) * t150 + (-t144 * t148 + t146 * t158) * t153, 0; 1, t149, -t147 * t154 * t146 + t149 * t148, 0, -t149 * t146 * t153 + (-t148 * t153 * t154 + t150 * t151) * t147, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:14:19
	% EndTime: 2019-10-10 12:14:19
	% DurationCPUTime: 0.08s
	% Computational Cost: add. (40->22), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->31)
	t198 = sin(pkin(7));
	t201 = cos(pkin(6));
	t217 = t198 * t201;
	t200 = cos(pkin(7));
	t206 = cos(qJ(2));
	t216 = t200 * t206;
	t199 = sin(pkin(6));
	t204 = sin(qJ(1));
	t215 = t204 * t199;
	t203 = sin(qJ(2));
	t214 = t204 * t203;
	t213 = t204 * t206;
	t207 = cos(qJ(1));
	t212 = t207 * t199;
	t211 = t207 * t203;
	t210 = t207 * t206;
	t191 = t201 * t210 - t214;
	t209 = -t191 * t200 + t198 * t212;
	t193 = -t201 * t213 - t211;
	t208 = t193 * t200 + t198 * t215;
	t205 = cos(qJ(3));
	t202 = sin(qJ(3));
	t197 = pkin(13) + qJ(5);
	t196 = cos(t197);
	t195 = sin(t197);
	t194 = -t201 * t214 + t210;
	t192 = t201 * t211 + t213;
	t190 = -t199 * t206 * t198 + t201 * t200;
	t189 = -t193 * t198 + t200 * t215;
	t188 = -t191 * t198 - t200 * t212;
	t1 = [0, t215, t189, 0, t194 * t202 - t208 * t205, (t194 * t205 + t208 * t202) * t195 - t189 * t196; 0, -t212, t188, 0, t192 * t202 + t209 * t205, (t192 * t205 - t209 * t202) * t195 - t188 * t196; 1, t201, t190, 0, -t205 * t217 + (t202 * t203 - t205 * t216) * t199, (t202 * t217 + (t202 * t216 + t203 * t205) * t199) * t195 - t190 * t196;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end