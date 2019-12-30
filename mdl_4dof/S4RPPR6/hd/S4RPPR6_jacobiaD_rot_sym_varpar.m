% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPPR6
% Use Code from Maple symbolic Code Generation
% 
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
% 
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% link_index [1x1 uint8]
%   Index des Segmentes, auf dem der Punkt C liegt.
%   Wie in S4RPPR6_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPR6_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR6_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR6_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR6_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR6_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:24
	% EndTime: 2019-12-29 12:45:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:25
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.04s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:24
	% EndTime: 2019-12-29 12:45:24
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:24
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.47s
	% Computational Cost: add. (213->24), mult. (456->71), div. (110->13), fcn. (596->7), ass. (0->33)
	t68 = sin(qJ(1));
	t80 = qJD(1) * t68;
	t67 = cos(pkin(6));
	t58 = 0.1e1 / t67 ^ 2;
	t66 = sin(pkin(6));
	t56 = t66 ^ 2;
	t61 = t68 ^ 2;
	t53 = t56 * t58 * t61 + 0.1e1;
	t69 = cos(qJ(1));
	t62 = t69 ^ 2;
	t84 = 0.1e1 / t53 ^ 2 * t62;
	t89 = t58 * t84;
	t49 = 0.1e1 / t53;
	t87 = (t49 - 0.1e1) * t66;
	t81 = t68 * t66;
	t48 = atan2(-t81, -t67);
	t46 = sin(t48);
	t47 = cos(t48);
	t44 = -t46 * t81 - t47 * t67;
	t41 = 0.1e1 / t44;
	t57 = 0.1e1 / t67;
	t42 = 0.1e1 / t44 ^ 2;
	t85 = t42 * t69;
	t83 = t56 * t57;
	t82 = t61 / t69 ^ 2;
	t79 = t57 * t89;
	t37 = (-t47 * t49 * t68 * t83 + t46 * t87) * t69;
	t55 = t66 * t56;
	t54 = t58 * t82 + 0.1e1;
	t43 = t41 * t42;
	t40 = t42 * t56 * t62 + 0.1e1;
	t36 = qJD(1) * t37;
	t1 = [(-t49 * t57 * t66 - 0.2e1 * t55 * t79) * t80, 0, 0, 0; (0.2e1 * (t37 * t85 + t41 * t68) / t40 ^ 2 * (-t36 * t43 * t62 - t80 * t85) * t56 + ((0.2e1 * t37 * t43 * t69 + t42 * t68) * t36 + (-t69 * t41 + ((t37 + (t55 * t89 + t87) * t69 * t46) * t68 - (0.2e1 * t61 * t56 ^ 2 * t79 + (t84 + (t61 - 0.2e1 * t62) * t49) * t83) * t69 * t47) * t42) * qJD(1)) / t40) * t66, 0, 0, 0; 0.2e1 * (-0.1e1 / t54 + (0.1e1 + t82) / t54 ^ 2 * t58) * (t61 / t62 + 0.1e1) / t69 * t57 * t80, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-12-29 12:45:25
	% EndTime: 2019-12-29 12:45:25
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end