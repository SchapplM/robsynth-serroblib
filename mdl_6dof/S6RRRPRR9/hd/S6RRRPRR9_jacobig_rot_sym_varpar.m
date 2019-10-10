% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRPRR9
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
% Datum: 2019-10-10 12:08
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRPRR9_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRPRR9_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
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
	% StartTime: 2019-10-10 12:08:07
	% EndTime: 2019-10-10 12:08:07
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (6->6), mult. (17->14), div. (0->0), fcn. (28->8), ass. (0->12)
	t105 = cos(pkin(6));
	t108 = cos(qJ(2));
	t112 = t105 * t108;
	t103 = sin(pkin(6));
	t107 = sin(qJ(1));
	t111 = t107 * t103;
	t109 = cos(qJ(1));
	t110 = t109 * t103;
	t106 = sin(qJ(2));
	t104 = cos(pkin(7));
	t102 = sin(pkin(7));
	t1 = [0, t111, -(-t109 * t106 - t107 * t112) * t102 + t104 * t111, 0, 0, 0; 0, -t110, -(-t107 * t106 + t109 * t112) * t102 - t104 * t110, 0, 0, 0; 1, t105, -t103 * t108 * t102 + t105 * t104, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:08
	% EndTime: 2019-10-10 12:08:08
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (25->16), mult. (72->33), div. (0->0), fcn. (105->12), ass. (0->25)
	t151 = sin(pkin(6));
	t157 = sin(qJ(1));
	t167 = t157 * t151;
	t156 = sin(qJ(2));
	t166 = t157 * t156;
	t159 = cos(qJ(2));
	t165 = t157 * t159;
	t160 = cos(qJ(1));
	t164 = t160 * t151;
	t163 = t160 * t156;
	t162 = t160 * t159;
	t149 = sin(pkin(13));
	t152 = cos(pkin(13));
	t155 = sin(qJ(3));
	t158 = cos(qJ(3));
	t161 = -t149 * t155 + t152 * t158;
	t154 = cos(pkin(6));
	t153 = cos(pkin(7));
	t150 = sin(pkin(7));
	t148 = -t158 * t149 - t155 * t152;
	t147 = -t154 * t165 - t163;
	t146 = t154 * t162 - t166;
	t145 = t161 * t153;
	t144 = t161 * t150;
	t1 = [0, t167, -t147 * t150 + t153 * t167, 0, -(-t154 * t166 + t162) * t148 - t147 * t145 - t144 * t167, 0; 0, -t164, -t146 * t150 - t153 * t164, 0, -(t154 * t163 + t165) * t148 - t146 * t145 + t144 * t164, 0; 1, t154, -t151 * t159 * t150 + t154 * t153, 0, -t154 * t144 + (-t145 * t159 - t148 * t156) * t151, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:08:09
	% EndTime: 2019-10-10 12:08:09
	% DurationCPUTime: 0.14s
	% Computational Cost: add. (52->25), mult. (148->51), div. (0->0), fcn. (211->14), ass. (0->34)
	t204 = sin(pkin(6));
	t211 = sin(qJ(1));
	t222 = t211 * t204;
	t210 = sin(qJ(2));
	t221 = t211 * t210;
	t214 = cos(qJ(2));
	t220 = t211 * t214;
	t215 = cos(qJ(1));
	t219 = t215 * t204;
	t218 = t215 * t210;
	t217 = t215 * t214;
	t202 = sin(pkin(13));
	t205 = cos(pkin(13));
	t209 = sin(qJ(3));
	t213 = cos(qJ(3));
	t216 = t213 * t202 + t209 * t205;
	t201 = -t209 * t202 + t213 * t205;
	t212 = cos(qJ(5));
	t208 = sin(qJ(5));
	t207 = cos(pkin(6));
	t206 = cos(pkin(7));
	t203 = sin(pkin(7));
	t199 = -t207 * t221 + t217;
	t198 = -t207 * t220 - t218;
	t197 = t207 * t218 + t220;
	t196 = t207 * t217 - t221;
	t195 = -t204 * t214 * t203 + t207 * t206;
	t194 = t216 * t206;
	t193 = t201 * t206;
	t192 = t216 * t203;
	t191 = t201 * t203;
	t190 = -t198 * t203 + t206 * t222;
	t189 = -t196 * t203 - t206 * t219;
	t1 = [0, t222, t190, 0, -t191 * t222 - t198 * t193 + t199 * t216, (t192 * t222 + t198 * t194 + t199 * t201) * t208 - t190 * t212; 0, -t219, t189, 0, t191 * t219 - t196 * t193 + t197 * t216, (-t192 * t219 + t196 * t194 + t197 * t201) * t208 - t189 * t212; 1, t207, t195, 0, -t207 * t191 + (-t193 * t214 + t210 * t216) * t204, (t207 * t192 + (t194 * t214 + t201 * t210) * t204) * t208 - t195 * t212;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end