% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RRRRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 12:54
% Revision: eb1f267a533306f0f157b6776e21de13647fd8af (2019-10-10)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRRRPR14_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RRRRPR14_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:12
	% EndTime: 2019-10-10 12:54:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (1->1), mult. (2->2), div. (0->0), fcn. (5->4), ass. (0->2)
	t43 = sin(pkin(6));
	t1 = [0, sin(qJ(1)) * t43, 0, 0, 0, 0; 0, -cos(qJ(1)) * t43, 0, 0, 0, 0; 1, cos(pkin(6)), 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.04s
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
	% StartTime: 2019-10-10 12:54:13
	% EndTime: 2019-10-10 12:54:13
	% DurationCPUTime: 0.05s
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
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (16->14), mult. (48->31), div. (0->0), fcn. (72->10), ass. (0->19)
	t179 = sin(pkin(6));
	t184 = sin(qJ(1));
	t193 = t184 * t179;
	t183 = sin(qJ(2));
	t192 = t184 * t183;
	t186 = cos(qJ(2));
	t191 = t184 * t186;
	t187 = cos(qJ(1));
	t190 = t187 * t179;
	t189 = t187 * t183;
	t188 = t187 * t186;
	t185 = cos(qJ(3));
	t182 = sin(qJ(3));
	t181 = cos(pkin(6));
	t180 = cos(pkin(7));
	t178 = sin(pkin(7));
	t177 = -t181 * t191 - t189;
	t176 = t181 * t188 - t192;
	t1 = [0, t193, -t177 * t178 + t180 * t193, (-t181 * t192 + t188) * t182 + (-t177 * t180 - t178 * t193) * t185, 0, 0; 0, -t190, -t176 * t178 - t180 * t190, (t181 * t189 + t191) * t182 + (-t176 * t180 + t178 * t190) * t185, 0, 0; 1, t181, -t179 * t186 * t178 + t181 * t180, -t181 * t178 * t185 + (-t180 * t185 * t186 + t182 * t183) * t179, 0, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 12:54:14
	% EndTime: 2019-10-10 12:54:14
	% DurationCPUTime: 0.12s
	% Computational Cost: add. (34->21), mult. (100->45), div. (0->0), fcn. (145->12), ass. (0->30)
	t189 = sin(pkin(7));
	t192 = cos(pkin(6));
	t210 = t189 * t192;
	t191 = cos(pkin(7));
	t199 = cos(qJ(2));
	t209 = t191 * t199;
	t190 = sin(pkin(6));
	t196 = sin(qJ(1));
	t208 = t196 * t190;
	t195 = sin(qJ(2));
	t207 = t196 * t195;
	t206 = t196 * t199;
	t200 = cos(qJ(1));
	t205 = t200 * t190;
	t204 = t200 * t195;
	t203 = t200 * t199;
	t185 = t192 * t203 - t207;
	t202 = -t185 * t191 + t189 * t205;
	t187 = -t192 * t206 - t204;
	t201 = t187 * t191 + t189 * t208;
	t198 = cos(qJ(3));
	t197 = cos(qJ(4));
	t194 = sin(qJ(3));
	t193 = sin(qJ(4));
	t188 = -t192 * t207 + t203;
	t186 = t192 * t204 + t206;
	t184 = -t189 * t190 * t199 + t191 * t192;
	t183 = -t187 * t189 + t191 * t208;
	t182 = -t185 * t189 - t191 * t205;
	t1 = [0, t208, t183, t188 * t194 - t198 * t201, 0, (t188 * t198 + t194 * t201) * t193 - t183 * t197; 0, -t205, t182, t186 * t194 + t198 * t202, 0, (t186 * t198 - t194 * t202) * t193 - t182 * t197; 1, t192, t184, -t198 * t210 + (t194 * t195 - t198 * t209) * t190, 0, (t194 * t210 + (t194 * t209 + t195 * t198) * t190) * t193 - t184 * t197;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end