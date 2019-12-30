% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPPPR3
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
% qJD [5x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S5RPPPR3_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 15:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPPPR3_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR3_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR3_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPPPR3_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR3_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (31->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:05
	% EndTime: 2019-12-29 15:46:06
	% DurationCPUTime: 0.46s
	% Computational Cost: add. (399->25), mult. (456->72), div. (110->13), fcn. (596->7), ass. (0->36)
	t71 = qJ(1) + pkin(7);
	t64 = sin(t71);
	t85 = qJD(1) * t64;
	t65 = cos(t71);
	t60 = t65 ^ 2;
	t72 = sin(pkin(8));
	t67 = t72 ^ 2;
	t93 = t60 * t67;
	t59 = t64 ^ 2;
	t73 = cos(pkin(8));
	t69 = 0.1e1 / t73 ^ 2;
	t88 = t59 * t69;
	t57 = t67 * t88 + 0.1e1;
	t55 = 0.1e1 / t57;
	t91 = (t55 - 0.1e1) * t72;
	t87 = t64 * t72;
	t54 = atan2(-t87, -t73);
	t49 = sin(t54);
	t50 = cos(t54);
	t47 = -t49 * t87 - t50 * t73;
	t44 = 0.1e1 / t47;
	t68 = 0.1e1 / t73;
	t45 = 0.1e1 / t47 ^ 2;
	t89 = t45 * t65;
	t86 = t67 * t68;
	t84 = t55 * t86;
	t56 = 0.1e1 / t57 ^ 2;
	t83 = t56 * t72 * t93;
	t40 = (-t50 * t64 * t84 + t49 * t91) * t65;
	t70 = t68 * t69;
	t62 = 0.1e1 / t65 ^ 2;
	t53 = t62 * t88 + 0.1e1;
	t46 = t44 * t45;
	t43 = t45 * t93 + 0.1e1;
	t39 = qJD(1) * t40;
	t1 = [(-t55 * t68 * t72 - 0.2e1 * t70 * t83) * t85, 0, 0, 0, 0; (0.2e1 * (t40 * t89 + t44 * t64) / t43 ^ 2 * (-t39 * t46 * t60 - t85 * t89) * t67 + ((0.2e1 * t40 * t46 * t65 + t45 * t64) * t39 + (-t65 * t44 + ((t40 + (t69 * t83 + t91) * t65 * t49) * t64 - (t59 * t84 + (-0.2e1 * t84 + (0.2e1 * t59 * t67 ^ 2 * t70 + t86) * t56) * t60) * t65 * t50) * t45) * qJD(1)) / t43) * t72, 0, 0, 0, 0; 0.2e1 * (-0.1e1 / t53 + (t59 * t62 + 0.1e1) / t53 ^ 2 * t69) * (t59 / t60 + 0.1e1) / t65 * t68 * t85, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 15:46:00
	% EndTime: 2019-12-29 15:46:00
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end