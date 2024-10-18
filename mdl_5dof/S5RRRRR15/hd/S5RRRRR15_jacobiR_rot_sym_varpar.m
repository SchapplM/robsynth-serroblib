% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
% 
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt. (0=Basis).
%   Siehe auch: S5RRRRR15_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5RRRRR15_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRRRR15_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (3->3), mult. (0->0), div. (0->0), fcn. (4->2), ass. (0->3)
	t9 = cos(qJ(1));
	t8 = sin(qJ(1));
	t1 = [-t8, 0, 0, 0, 0; t9, 0, 0, 0, 0; 0, 0, 0, 0, 0; -t9, 0, 0, 0, 0; -t8, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (9->7), mult. (28->12), div. (0->0), fcn. (48->6), ass. (0->15)
	t50 = sin(qJ(2));
	t51 = sin(qJ(1));
	t57 = t51 * t50;
	t52 = cos(qJ(2));
	t56 = t51 * t52;
	t53 = cos(qJ(1));
	t55 = t53 * t50;
	t54 = t53 * t52;
	t49 = cos(pkin(5));
	t48 = sin(pkin(5));
	t47 = -t49 * t57 + t54;
	t46 = -t49 * t56 - t55;
	t45 = -t49 * t55 - t56;
	t44 = -t49 * t54 + t57;
	t1 = [t45, t46, 0, 0, 0; t47, -t44, 0, 0, 0; 0, t48 * t52, 0, 0, 0; t44, -t47, 0, 0, 0; t46, t45, 0, 0, 0; 0, -t48 * t50, 0, 0, 0; t53 * t48, 0, 0, 0, 0; t51 * t48, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (104->15), mult. (46->16), div. (0->0), fcn. (72->9), ass. (0->23)
	t92 = sin(qJ(1));
	t97 = -t92 / 0.2e1;
	t93 = cos(qJ(1));
	t96 = -t93 / 0.2e1;
	t90 = qJ(2) + qJ(3);
	t86 = pkin(5) + t90;
	t84 = cos(t86);
	t87 = pkin(5) - t90;
	t85 = cos(t87);
	t95 = t85 + t84;
	t82 = sin(t86);
	t83 = sin(t87);
	t81 = t82 - t83;
	t89 = cos(t90);
	t75 = t81 * t96 - t92 * t89;
	t94 = t81 * t97 + t93 * t89;
	t91 = sin(pkin(5));
	t88 = sin(t90);
	t80 = t84 / 0.2e1 - t85 / 0.2e1;
	t79 = t82 / 0.2e1 + t83 / 0.2e1;
	t76 = -t93 * t88 + t95 * t97;
	t74 = t92 * t88 + t95 * t96;
	t1 = [t75, t76, t76, 0, 0; t94, -t74, -t74, 0, 0; 0, t79, t79, 0, 0; t74, -t94, -t94, 0, 0; t76, t75, t75, 0, 0; 0, t80, t80, 0, 0; t93 * t91, 0, 0, 0, 0; t92 * t91, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:07
	% EndTime: 2024-09-27 22:28:07
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (202->15), mult. (78->16), div. (0->0), fcn. (96->9), ass. (0->21)
	t96 = qJ(2) + qJ(3) + qJ(4);
	t93 = pkin(5) + t96;
	t105 = cos(t93) / 0.2e1;
	t104 = sin(t93) / 0.2e1;
	t103 = pkin(5) - t96;
	t101 = sin(t103);
	t89 = t104 - t101 / 0.2e1;
	t95 = cos(t96);
	t98 = sin(qJ(1));
	t99 = cos(qJ(1));
	t84 = -t89 * t99 - t95 * t98;
	t86 = t89 * t98 - t95 * t99;
	t102 = cos(t103);
	t100 = t102 / 0.2e1 + t105;
	t97 = sin(pkin(5));
	t94 = sin(t96);
	t90 = t105 - t102 / 0.2e1;
	t88 = t104 + t101 / 0.2e1;
	t85 = -t100 * t98 - t99 * t94;
	t83 = -t100 * t99 + t98 * t94;
	t1 = [t84, t85, t85, t85, 0; -t86, -t83, -t83, -t83, 0; 0, t88, t88, t88, 0; t83, t86, t86, t86, 0; t85, t84, t84, t84, 0; 0, t90, t90, t90, 0; t99 * t97, 0, 0, 0, 0; t98 * t97, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-27 22:28:08
	% EndTime: 2024-09-27 22:28:08
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (238->35), mult. (318->74), div. (0->0), fcn. (445->10), ass. (0->50)
	t185 = sin(pkin(6));
	t188 = cos(pkin(5));
	t211 = t185 * t188;
	t186 = sin(pkin(5));
	t190 = sin(qJ(1));
	t210 = t186 * t190;
	t192 = cos(qJ(1));
	t209 = t186 * t192;
	t187 = cos(pkin(6));
	t189 = sin(qJ(5));
	t208 = t187 * t189;
	t191 = cos(qJ(5));
	t207 = t187 * t191;
	t206 = t188 * t190;
	t205 = t188 * t192;
	t204 = t190 * t189;
	t203 = t191 * t190;
	t202 = t192 * t189;
	t201 = t192 * t191;
	t200 = t185 * t210;
	t199 = t185 * t209;
	t198 = t188 * t202;
	t197 = t188 * t204;
	t196 = t188 * t203;
	t195 = t188 * t201;
	t172 = t187 * t197 - t201;
	t178 = t187 * t202 + t196;
	t184 = qJ(2) + qJ(3) + qJ(4);
	t182 = sin(t184);
	t183 = cos(t184);
	t194 = -t172 * t183 - t178 * t182 + t189 * t200;
	t175 = t187 * t195 - t204;
	t177 = t187 * t203 + t198;
	t193 = -t175 * t183 + t177 * t182 + t191 * t199;
	t180 = t186 * t182 * t185;
	t179 = -t187 * t201 + t197;
	t176 = t187 * t204 - t195;
	t174 = t187 * t198 + t203;
	t173 = -t187 * t196 - t202;
	t171 = (-t182 * t208 + t183 * t191) * t186;
	t170 = (-t182 * t207 - t183 * t189) * t186;
	t169 = (t182 * t205 + t183 * t190) * t185;
	t168 = (-t182 * t206 + t183 * t192) * t185;
	t167 = -t175 * t182 - t177 * t183;
	t166 = -t174 * t182 - t176 * t183;
	t165 = -t173 * t182 + t179 * t183;
	t164 = t172 * t182 - t178 * t183;
	t163 = -t174 * t183 + t176 * t182 + t189 * t199;
	t162 = t173 * t183 + t179 * t182 + t191 * t200;
	t1 = [t163, t164, t164, t164, t162; t194, t166, t166, t166, -t193; 0, t171, t171, t171, t191 * t211 + (-t182 * t189 + t183 * t207) * t186; t193, t165, t165, t165, -t194; t162, t167, t167, t167, t163; 0, t170, t170, t170, -t189 * t211 + (-t182 * t191 - t183 * t208) * t186; t187 * t209 + (-t182 * t190 + t183 * t205) * t185, t168, t168, t168, 0; t187 * t210 + (t182 * t192 + t183 * t206) * t185, t169, t169, t169, 0; 0, t180, t180, t180, 0;];
	JR_rot = t1;
end