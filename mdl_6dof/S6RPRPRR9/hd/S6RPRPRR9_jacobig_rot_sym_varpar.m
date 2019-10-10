% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für beliebiges Segment von
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
% 
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-10 01:00
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RPRPRR9_jacobig_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),uint8(0),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobig_rot_sym_varpar: qJ has to be [6x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S6RPRPRR9_jacobig_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobig_rot_sym_varpar: pkin has to be [13x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobig_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobig_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:09
	% EndTime: 2019-10-10 01:00:09
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobig_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 1, 0, 0, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobig_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t78 = sin(pkin(6));
	t80 = cos(pkin(7));
	t85 = t78 * t80;
	t79 = cos(pkin(12));
	t81 = cos(pkin(6));
	t84 = t79 * t81;
	t83 = cos(qJ(1));
	t82 = sin(qJ(1));
	t77 = sin(pkin(7));
	t76 = sin(pkin(12));
	t1 = [0, 0, -(-t83 * t76 - t82 * t84) * t77 + t82 * t85, 0, 0, 0; 0, 0, -(-t82 * t76 + t83 * t84) * t77 - t83 * t85, 0, 0, 0; 1, 0, -t78 * t79 * t77 + t81 * t80, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobig_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (5->5), mult. (15->13), div. (0->0), fcn. (23->8), ass. (0->11)
	t95 = sin(pkin(6));
	t97 = cos(pkin(7));
	t102 = t95 * t97;
	t96 = cos(pkin(12));
	t98 = cos(pkin(6));
	t101 = t96 * t98;
	t100 = cos(qJ(1));
	t99 = sin(qJ(1));
	t94 = sin(pkin(7));
	t93 = sin(pkin(12));
	t1 = [0, 0, -(-t100 * t93 - t99 * t101) * t94 + t99 * t102, 0, 0, 0; 0, 0, -(t100 * t101 - t99 * t93) * t94 - t100 * t102, 0, 0, 0; 1, 0, -t95 * t96 * t94 + t98 * t97, 0, 0, 0;];
	Jg_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobig_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:10
	% EndTime: 2019-10-10 01:00:10
	% DurationCPUTime: 0.09s
	% Computational Cost: add. (24->15), mult. (70->33), div. (0->0), fcn. (100->12), ass. (0->25)
	t137 = sin(pkin(6));
	t143 = sin(qJ(1));
	t152 = t137 * t143;
	t145 = cos(qJ(1));
	t151 = t137 * t145;
	t135 = sin(pkin(12));
	t150 = t143 * t135;
	t139 = cos(pkin(12));
	t149 = t143 * t139;
	t148 = t145 * t135;
	t147 = t145 * t139;
	t134 = sin(pkin(13));
	t138 = cos(pkin(13));
	t142 = sin(qJ(3));
	t144 = cos(qJ(3));
	t146 = -t134 * t142 + t138 * t144;
	t141 = cos(pkin(6));
	t140 = cos(pkin(7));
	t136 = sin(pkin(7));
	t133 = -t134 * t144 - t138 * t142;
	t132 = -t141 * t149 - t148;
	t131 = t141 * t147 - t150;
	t130 = t146 * t140;
	t129 = t146 * t136;
	t1 = [0, 0, -t132 * t136 + t140 * t152, 0, -(-t141 * t150 + t147) * t133 - t132 * t130 - t129 * t152, 0; 0, 0, -t131 * t136 - t140 * t151, 0, -(t141 * t148 + t149) * t133 - t131 * t130 + t129 * t151, 0; 1, 0, -t136 * t137 * t139 + t140 * t141, 0, -t141 * t129 + (-t130 * t139 - t133 * t135) * t137, 0;];
	Jg_rot = t1;
elseif link_index == 6
	%% Symbolic Calculation
	% From jacobig_rot_6_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-10 01:00:11
	% EndTime: 2019-10-10 01:00:11
	% DurationCPUTime: 0.11s
	% Computational Cost: add. (51->24), mult. (146->51), div. (0->0), fcn. (206->14), ass. (0->34)
	t182 = sin(pkin(6));
	t189 = sin(qJ(1));
	t199 = t182 * t189;
	t192 = cos(qJ(1));
	t198 = t182 * t192;
	t180 = sin(pkin(12));
	t197 = t189 * t180;
	t184 = cos(pkin(12));
	t196 = t189 * t184;
	t195 = t192 * t180;
	t194 = t192 * t184;
	t179 = sin(pkin(13));
	t183 = cos(pkin(13));
	t188 = sin(qJ(3));
	t191 = cos(qJ(3));
	t193 = t191 * t179 + t188 * t183;
	t178 = -t188 * t179 + t191 * t183;
	t190 = cos(qJ(5));
	t187 = sin(qJ(5));
	t186 = cos(pkin(6));
	t185 = cos(pkin(7));
	t181 = sin(pkin(7));
	t176 = -t186 * t197 + t194;
	t175 = -t186 * t196 - t195;
	t174 = t186 * t195 + t196;
	t173 = t186 * t194 - t197;
	t172 = -t182 * t184 * t181 + t186 * t185;
	t171 = t193 * t185;
	t170 = t178 * t185;
	t169 = t193 * t181;
	t168 = t178 * t181;
	t167 = -t175 * t181 + t185 * t199;
	t166 = -t173 * t181 - t185 * t198;
	t1 = [0, 0, t167, 0, -t168 * t199 - t175 * t170 + t176 * t193, (t169 * t199 + t175 * t171 + t176 * t178) * t187 - t167 * t190; 0, 0, t166, 0, t168 * t198 - t173 * t170 + t174 * t193, (-t169 * t198 + t173 * t171 + t174 * t178) * t187 - t166 * t190; 1, 0, t172, 0, -t186 * t168 + (-t170 * t184 + t180 * t193) * t182, (t186 * t169 + (t171 * t184 + t178 * t180) * t182) * t187 - t172 * t190;];
	Jg_rot = t1;
else
	Jg_rot=NaN(3,6);
end