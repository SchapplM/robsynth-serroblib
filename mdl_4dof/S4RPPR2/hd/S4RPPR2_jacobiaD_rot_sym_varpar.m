% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4RPPR2
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
%   Wie in S4RPPR2_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:41
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPPR2_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4RPPR2_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.02s
	% Computational Cost: add. (9->0), mult. (37->0), div. (15->0), fcn. (22->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (52->9), mult. (138->18), div. (22->4), fcn. (156->4), ass. (0->15)
	t42 = sin(pkin(6));
	t43 = cos(pkin(6));
	t44 = sin(qJ(1));
	t45 = cos(qJ(1));
	t37 = -t44 * t42 - t45 * t43;
	t35 = 0.1e1 / t37 ^ 2;
	t38 = t45 * t42 - t44 * t43;
	t50 = t38 ^ 2 * t35;
	t51 = t35 * t37;
	t32 = t38 * qJD(1);
	t34 = 0.1e1 / t37;
	t49 = t34 * t50;
	t48 = t32 * t51;
	t30 = 0.1e1 + t50;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t34 * t37 + t50) / t30 ^ 2 * (t32 * t49 + t48) + (-0.2e1 * t48 + (t34 - 0.2e1 * t49 - t51) * t32) / t30, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:41:54
	% EndTime: 2019-10-09 20:41:54
	% DurationCPUTime: 0.13s
	% Computational Cost: add. (306->14), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->22)
	t70 = pkin(6) + qJ(4);
	t60 = sin(t70);
	t61 = cos(qJ(1));
	t67 = cos(t70);
	t75 = sin(qJ(1));
	t65 = t75 * t60 + t61 * t67;
	t50 = 0.1e1 / t65 ^ 2;
	t77 = t50 * t65;
	t76 = qJD(1) - qJD(4);
	t49 = 0.1e1 / t65;
	t64 = -t61 * t60 + t75 * t67;
	t44 = t76 * t64;
	t48 = t64 ^ 2;
	t72 = t48 * t50;
	t47 = 0.1e1 + t72;
	t73 = t76 * t77;
	t69 = t64 * t73;
	t51 = t49 * t50;
	t71 = t48 * t51;
	t74 = (t44 * t71 + t69) / t47 ^ 2;
	t45 = 0.1e1 / t47;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0.2e1 * (t49 * t65 + t72) * t74 + (-0.2e1 * t69 - (-t49 + 0.2e1 * t71 + t77) * t44) * t45, 0, 0, -0.2e1 * t74 - 0.2e1 * (-t45 * t73 - (t44 * t45 * t51 - t50 * t74) * t64) * t64;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end