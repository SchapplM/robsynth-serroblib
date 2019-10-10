% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RRPR2
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
%   Wie in S4RRPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:51
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RRPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RRPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_jacobiaD_rot_sym_varpar: pkin has to be [5x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:50
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.03s
	% Computational Cost: add. (74->0), mult. (74->0), div. (30->0), fcn. (44->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:51:51
	% EndTime: 2019-10-09 20:51:51
	% DurationCPUTime: 0.15s
	% Computational Cost: add. (514->14), mult. (492->30), div. (62->4), fcn. (544->4), ass. (0->23)
	t79 = qJ(1) + qJ(2);
	t68 = cos(t79);
	t70 = sin(qJ(4));
	t76 = sin(t79);
	t84 = cos(qJ(4));
	t74 = t68 * t84 + t76 * t70;
	t58 = 0.1e1 / t74 ^ 2;
	t86 = t58 * t74;
	t85 = qJD(4) - qJD(1) - qJD(2);
	t57 = 0.1e1 / t74;
	t73 = -t68 * t70 + t76 * t84;
	t52 = t85 * t73;
	t56 = t73 ^ 2;
	t81 = t56 * t58;
	t55 = 0.1e1 + t81;
	t82 = t85 * t86;
	t78 = t73 * t82;
	t59 = t57 * t58;
	t80 = t56 * t59;
	t83 = (-t52 * t80 - t78) / t55 ^ 2;
	t53 = 0.1e1 / t55;
	t49 = 0.2e1 * (t57 * t74 + t81) * t83 + (0.2e1 * t78 + (-t57 + 0.2e1 * t80 + t86) * t52) * t53;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; t49, t49, 0, -0.2e1 * t83 - 0.2e1 * (t53 * t82 - (-t52 * t53 * t59 - t58 * t83) * t73) * t73;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end