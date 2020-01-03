% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RRPRP5_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5RRPRP5_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RRPRP5_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (3->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 1, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 19:55:55
	% EndTime: 2019-12-31 19:55:55
	% DurationCPUTime: 0.07s
	% Computational Cost: add. (471->18), mult. (227->41), div. (75->9), fcn. (351->7), ass. (0->29)
	t57 = cos(qJ(1));
	t58 = t57 ^ 2;
	t52 = qJ(2) + pkin(8) + qJ(4);
	t51 = cos(t52);
	t50 = sin(t52);
	t56 = sin(qJ(1));
	t60 = t56 * t50;
	t44 = atan2(-t60, -t51);
	t42 = sin(t44);
	t43 = cos(t44);
	t40 = -t42 * t60 - t43 * t51;
	t39 = 0.1e1 / t40 ^ 2;
	t65 = t39 * t50;
	t64 = t42 * t51;
	t47 = t50 ^ 2;
	t49 = 0.1e1 / t51 ^ 2;
	t63 = t47 * t49;
	t53 = t56 ^ 2;
	t62 = t53 / t58;
	t45 = 0.1e1 / (t53 * t63 + 0.1e1);
	t61 = t56 * t45;
	t46 = 0.1e1 / (t49 * t62 + 0.1e1);
	t59 = 0.1e1 / t57 * t49 * t46 * t60;
	t48 = 0.1e1 / t51;
	t41 = (0.1e1 + t63) * t61;
	t38 = 0.1e1 / t40;
	t37 = 0.1e1 / (t58 * t47 * t39 + 0.1e1);
	t36 = (t51 * t38 - (-t56 * t64 + t43 * t50 + (-t43 * t60 + t64) * t41) * t65) * t57 * t37;
	t1 = [t57 * t50 * t48 * t45, t41, 0, t41, 0; (-t38 * t60 - (-t43 * t47 * t48 * t61 + (t45 - 0.1e1) * t50 * t42) * t58 * t65) * t37, t36, 0, t36, 0; (-0.1e1 - t62) * t48 * t46, -t59, 0, -t59, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end