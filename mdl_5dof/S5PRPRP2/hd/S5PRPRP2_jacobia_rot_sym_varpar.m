% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5PRPRP2
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
%   Wie in S5PRPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% Ja_rot [3x5]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja_rot = S5PRPRP2_jacobia_rot_sym_varpar(qJ, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_jacobia_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5PRPRP2_jacobia_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_jacobia_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobia_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobia_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobia_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (9->0), mult. (6->0), div. (5->0), fcn. (6->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 1, 0, 0, 0;];
	Ja_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobia_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN; NaN, NaN, NaN, NaN, NaN;];
	Ja_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobia_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (115->15), mult. (116->36), div. (25->9), fcn. (175->9), ass. (0->24)
	t36 = cos(pkin(8));
	t34 = pkin(7) + qJ(2);
	t30 = sin(t34);
	t35 = sin(pkin(8));
	t42 = t30 * t35;
	t28 = atan2(-t42, -t36);
	t26 = sin(t28);
	t27 = cos(t28);
	t21 = -t26 * t42 - t27 * t36;
	t31 = cos(t34);
	t44 = 0.1e1 / t21 ^ 2 * t31 ^ 2;
	t32 = t35 ^ 2;
	t29 = 0.1e1 / (0.1e1 + t30 ^ 2 * t32 / t36 ^ 2);
	t43 = t29 / t36;
	t37 = sin(qJ(4));
	t41 = t36 * t37;
	t38 = cos(qJ(4));
	t40 = t36 * t38;
	t25 = t30 * t37 + t31 * t40;
	t23 = 0.1e1 / t25 ^ 2;
	t24 = -t30 * t38 + t31 * t41;
	t39 = t24 ^ 2 * t23 + 0.1e1;
	t22 = 0.1e1 / t39;
	t1 = [0, t31 * t35 * t43, 0, 0, 0; 0, (-0.1e1 / t21 * t42 - (-t27 * t30 * t32 * t43 + (t29 - 0.1e1) * t35 * t26) * t35 * t44) / (t32 * t44 + 0.1e1), 0, 0, 0; 0, ((-t30 * t41 - t31 * t38) / t25 - (-t30 * t40 + t31 * t37) * t24 * t23) * t22, 0, t39 * t22, 0;];
	Ja_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobia_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-05 15:31:30
	% EndTime: 2019-12-05 15:31:30
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (115->15), mult. (116->36), div. (25->9), fcn. (175->9), ass. (0->24)
	t38 = cos(pkin(8));
	t36 = pkin(7) + qJ(2);
	t32 = sin(t36);
	t37 = sin(pkin(8));
	t44 = t32 * t37;
	t30 = atan2(-t44, -t38);
	t28 = sin(t30);
	t29 = cos(t30);
	t23 = -t28 * t44 - t29 * t38;
	t33 = cos(t36);
	t46 = 0.1e1 / t23 ^ 2 * t33 ^ 2;
	t34 = t37 ^ 2;
	t31 = 0.1e1 / (0.1e1 + t32 ^ 2 * t34 / t38 ^ 2);
	t45 = t31 / t38;
	t39 = sin(qJ(4));
	t43 = t38 * t39;
	t40 = cos(qJ(4));
	t42 = t38 * t40;
	t27 = t32 * t39 + t33 * t42;
	t25 = 0.1e1 / t27 ^ 2;
	t26 = -t32 * t40 + t33 * t43;
	t41 = t26 ^ 2 * t25 + 0.1e1;
	t24 = 0.1e1 / t41;
	t1 = [0, t33 * t37 * t45, 0, 0, 0; 0, (-0.1e1 / t23 * t44 - (-t29 * t32 * t34 * t45 + (t31 - 0.1e1) * t37 * t28) * t37 * t46) / (t34 * t46 + 0.1e1), 0, 0, 0; 0, ((-t32 * t43 - t33 * t40) / t27 - (-t32 * t42 + t33 * t39) * t26 * t25) * t24, 0, t41 * t24, 0;];
	Ja_rot = t1;
else
	Ja_rot=NaN(3,5);
end