% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S5RPRPR10
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
%   Wie in S5RPRPR10_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% JaD_rot [3x5]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S5RPRPR10_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),uint8(0),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_jacobiaD_rot_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_jacobiaD_rot_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S5RPRPR10_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_jacobiaD_rot_sym_varpar: pkin has to be [8x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (128->13), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->21)
	t54 = sin(qJ(3));
	t55 = cos(qJ(1));
	t67 = sin(qJ(1));
	t68 = cos(qJ(3));
	t59 = t67 * t54 + t55 * t68;
	t44 = 0.1e1 / t59 ^ 2;
	t70 = t44 * t59;
	t69 = qJD(1) - qJD(3);
	t43 = 0.1e1 / t59;
	t58 = -t54 * t55 + t67 * t68;
	t38 = t69 * t58;
	t42 = t58 ^ 2;
	t64 = t42 * t44;
	t41 = 0.1e1 + t64;
	t65 = t69 * t70;
	t62 = t58 * t65;
	t45 = t43 * t44;
	t63 = t42 * t45;
	t66 = (t38 * t63 + t62) / t41 ^ 2;
	t39 = 0.1e1 / t41;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0.2e1 * (t43 * t59 + t64) * t66 + (-0.2e1 * t62 - (-t43 + 0.2e1 * t63 + t70) * t38) * t39, 0, -0.2e1 * t66 - 0.2e1 * (-t39 * t65 - (t38 * t39 * t45 - t44 * t66) * t58) * t58, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (306->14), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->22)
	t71 = qJ(3) + pkin(8);
	t61 = sin(t71);
	t62 = cos(qJ(1));
	t68 = cos(t71);
	t76 = sin(qJ(1));
	t66 = t76 * t61 + t62 * t68;
	t51 = 0.1e1 / t66 ^ 2;
	t78 = t51 * t66;
	t77 = qJD(1) - qJD(3);
	t50 = 0.1e1 / t66;
	t65 = -t62 * t61 + t76 * t68;
	t45 = t77 * t65;
	t49 = t65 ^ 2;
	t73 = t49 * t51;
	t48 = 0.1e1 + t73;
	t74 = t77 * t78;
	t70 = t65 * t74;
	t52 = t50 * t51;
	t72 = t49 * t52;
	t75 = (t45 * t72 + t70) / t48 ^ 2;
	t46 = 0.1e1 / t48;
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0.2e1 * (t50 * t66 + t73) * t75 + (-0.2e1 * t70 - (-t50 + 0.2e1 * t72 + t78) * t45) * t46, 0, -0.2e1 * t75 - 0.2e1 * (-t46 * t74 - (t45 * t46 * t52 - t51 * t75) * t65) * t65, 0, 0;];
	JaD_rot = t1;
elseif link_index == 5
	%% Symbolic Calculation
	% From jacobiaD_rot_5_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-31 18:26:21
	% EndTime: 2019-12-31 18:26:21
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0, 0; 0, 0, 0, 0, 0; 0, 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,5);
end