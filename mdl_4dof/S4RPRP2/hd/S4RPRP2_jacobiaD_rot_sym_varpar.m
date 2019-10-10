% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPRP2
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
%   Wie in S4RPRP2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:44
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPRP2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPRP2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_jacobiaD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:44:44
	% EndTime: 2019-10-09 20:44:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:44:44
	% EndTime: 2019-10-09 20:44:44
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:44:44
	% EndTime: 2019-10-09 20:44:44
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:44:44
	% EndTime: 2019-10-09 20:44:44
	% DurationCPUTime: 0.15s
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
	t58 = -t55 * t54 + t67 * t68;
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
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t43 * t59 + t64) * t66 + (-0.2e1 * t62 - (-t43 + 0.2e1 * t63 + t70) * t38) * t39, 0, -0.2e1 * t66 - 0.2e1 * (-t39 * t65 - (t38 * t39 * t45 - t44 * t66) * t58) * t58, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:44:44
	% EndTime: 2019-10-09 20:44:44
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (128->13), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->21)
	t56 = sin(qJ(3));
	t57 = cos(qJ(1));
	t69 = sin(qJ(1));
	t70 = cos(qJ(3));
	t61 = t69 * t56 + t57 * t70;
	t46 = 0.1e1 / t61 ^ 2;
	t72 = t46 * t61;
	t71 = qJD(1) - qJD(3);
	t45 = 0.1e1 / t61;
	t60 = -t56 * t57 + t69 * t70;
	t40 = t71 * t60;
	t44 = t60 ^ 2;
	t66 = t44 * t46;
	t43 = 0.1e1 + t66;
	t67 = t71 * t72;
	t64 = t60 * t67;
	t47 = t45 * t46;
	t65 = t44 * t47;
	t68 = (t40 * t65 + t64) / t43 ^ 2;
	t41 = 0.1e1 / t43;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t45 * t61 + t66) * t68 + (-0.2e1 * t64 - (-t45 + 0.2e1 * t65 + t72) * t40) * t41, 0, -0.2e1 * t68 - 0.2e1 * (-t41 * t67 - (t40 * t41 * t47 - t46 * t68) * t60) * t60, 0;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end