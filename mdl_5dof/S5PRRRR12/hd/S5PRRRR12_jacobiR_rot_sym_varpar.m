% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für beliebiges Segment von
% S5PRRRR12
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
%   Siehe auch: S5PRRRR12_fkine_fixb_rotmat_mdh_sym_varpar.m
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% JR_rot [9x5]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S5PRRRR12_jacobiR_rot_sym_varpar(qJ, link_index, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_jacobiR_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRRRR12_jacobiR_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_jacobiR_rot_sym_varpar: pkin has to be [11x1] (double)');
JR_rot=NaN(9,5);
if link_index == 0
	%% Symbolic Calculation
	% From jacobiR_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiR_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiR_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (5->5), mult. (14->12), div. (0->0), fcn. (24->6), ass. (0->9)
	t22 = cos(pkin(5));
	t23 = sin(qJ(2));
	t26 = t22 * t23;
	t24 = cos(qJ(2));
	t25 = t22 * t24;
	t21 = cos(pkin(11));
	t20 = sin(pkin(5));
	t19 = sin(pkin(11));
	t1 = [0, -t19 * t25 - t21 * t23, 0, 0, 0; 0, -t19 * t23 + t21 * t25, 0, 0, 0; 0, t20 * t24, 0, 0, 0; 0, t19 * t26 - t21 * t24, 0, 0, 0; 0, -t19 * t24 - t21 * t26, 0, 0, 0; 0, -t20 * t23, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiR_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (76->11), mult. (40->14), div. (0->0), fcn. (48->8), ass. (0->20)
	t60 = qJ(2) + qJ(3);
	t62 = cos(pkin(11));
	t61 = sin(pkin(11));
	t59 = cos(t60);
	t58 = sin(t60);
	t57 = pkin(5) - t60;
	t56 = pkin(5) + t60;
	t55 = cos(t57);
	t54 = sin(t56);
	t53 = cos(t56) / 0.2e1;
	t52 = sin(t57) / 0.2e1;
	t51 = t55 / 0.2e1 + t53;
	t50 = t53 - t55 / 0.2e1;
	t49 = t52 - t54 / 0.2e1;
	t48 = t54 / 0.2e1 + t52;
	t47 = -t61 * t49 - t62 * t59;
	t46 = -t61 * t51 - t62 * t58;
	t45 = t62 * t49 - t61 * t59;
	t44 = t62 * t51 - t61 * t58;
	t1 = [0, t46, t46, 0, 0; 0, t44, t44, 0, 0; 0, t48, t48, 0, 0; 0, t47, t47, 0, 0; 0, t45, t45, 0, 0; 0, t50, t50, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiR_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:12
	% DurationCPUTime: 0.00s
	% Computational Cost: add. (162->11), mult. (60->14), div. (0->0), fcn. (72->8), ass. (0->20)
	t72 = qJ(2) + qJ(3) + qJ(4);
	t74 = cos(pkin(11));
	t73 = sin(pkin(11));
	t71 = cos(t72);
	t70 = sin(t72);
	t69 = pkin(5) - t72;
	t68 = pkin(5) + t72;
	t67 = cos(t69);
	t66 = sin(t68);
	t65 = cos(t68) / 0.2e1;
	t64 = sin(t69) / 0.2e1;
	t63 = t67 / 0.2e1 + t65;
	t62 = t65 - t67 / 0.2e1;
	t61 = t64 - t66 / 0.2e1;
	t60 = t66 / 0.2e1 + t64;
	t59 = -t73 * t61 - t74 * t71;
	t58 = -t73 * t63 - t74 * t70;
	t57 = t74 * t61 - t73 * t71;
	t56 = t74 * t63 - t73 * t70;
	t1 = [0, t58, t58, t58, 0; 0, t56, t56, t56, 0; 0, t60, t60, t60, 0; 0, t59, t59, t59, 0; 0, t57, t57, t57, 0; 0, t62, t62, t62, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JR_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiR_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2024-09-28 18:09:12
	% EndTime: 2024-09-28 18:09:13
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (194->29), mult. (252->67), div. (0->0), fcn. (353->10), ass. (0->39)
	t178 = qJ(2) + qJ(3) + qJ(4);
	t176 = sin(t178);
	t184 = cos(pkin(5));
	t196 = t176 * t184;
	t180 = sin(pkin(6));
	t181 = sin(pkin(5));
	t195 = t180 * t181;
	t183 = cos(pkin(6));
	t185 = sin(qJ(5));
	t194 = t183 * t185;
	t186 = cos(qJ(5));
	t193 = t183 * t186;
	t192 = t184 * t185;
	t191 = t184 * t186;
	t190 = t185 * t195;
	t189 = t186 * t195;
	t188 = t183 * t192;
	t187 = t183 * t191;
	t182 = cos(pkin(11));
	t179 = sin(pkin(11));
	t177 = cos(t178);
	t175 = t176 * t195;
	t174 = t179 * t192 - t182 * t193;
	t173 = t179 * t191 + t182 * t194;
	t172 = -t179 * t193 - t182 * t192;
	t171 = t179 * t194 - t182 * t191;
	t170 = -t179 * t185 + t182 * t187;
	t169 = t179 * t186 + t182 * t188;
	t168 = -t179 * t187 - t182 * t185;
	t167 = t179 * t188 - t182 * t186;
	t166 = (-t176 * t194 + t177 * t186) * t181;
	t165 = (-t176 * t193 - t177 * t185) * t181;
	t164 = (t177 * t179 + t182 * t196) * t180;
	t163 = (t177 * t182 - t179 * t196) * t180;
	t162 = -t170 * t176 + t172 * t177;
	t161 = -t169 * t176 - t171 * t177;
	t160 = -t168 * t176 + t174 * t177;
	t159 = t167 * t176 - t173 * t177;
	t1 = [0, t159, t159, t159, t168 * t177 + t174 * t176 + t179 * t189; 0, t161, t161, t161, t170 * t177 + t172 * t176 - t182 * t189; 0, t166, t166, t166, t180 * t191 + (-t176 * t185 + t177 * t193) * t181; 0, t160, t160, t160, t167 * t177 + t173 * t176 - t179 * t190; 0, t162, t162, t162, -t169 * t177 + t171 * t176 + t182 * t190; 0, t165, t165, t165, -t180 * t192 + (-t176 * t186 - t177 * t194) * t181; 0, t163, t163, t163, 0; 0, t164, t164, t164, 0; 0, t175, t175, t175, 0;];
	JR_rot = t1;
end