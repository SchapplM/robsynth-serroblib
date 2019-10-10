% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für beliebiges Segment von
% S4PPRR1
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
%   Wie in S4PPRR1_fkine_fixb_rotmat_mdh_sym_varpar.m (1=Basis).
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
% 
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-10-09 20:23
% Revision: ee6bc4d0f60ba4b3bab3f447780ef990a2753b00 (2019-10-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPRR1_jacobiaD_rot_sym_varpar(qJ, qJD, link_index, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),uint8(0),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_jacobiaD_rot_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_jacobiaD_rot_sym_varpar: qJD has to be [4x1] (double)');
assert(isa(link_index,'uint8') && all(size(link_index) == [1 1]), ...
	'S4PPRR1_jacobiaD_rot_sym_varpar: link_index has to be [1x1] uint8');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_jacobiaD_rot_sym_varpar: pkin has to be [6x1] (double)');
if link_index == 0
	%% Symbolic Calculation
	% From jacobiaD_rot_0_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 1
	%% Symbolic Calculation
	% From jacobiaD_rot_1_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 2
	%% Symbolic Calculation
	% From jacobiaD_rot_2_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.01s
	% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, 0;];
	JaD_rot = t1;
elseif link_index == 3
	%% Symbolic Calculation
	% From jacobiaD_rot_3_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.06s
	% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
	t33 = sin(pkin(6));
	t34 = cos(pkin(6));
	t35 = sin(qJ(3));
	t36 = cos(qJ(3));
	t31 = t33 * t35 + t34 * t36;
	t28 = 0.1e1 / t31 ^ 2;
	t43 = qJD(3) * t28;
	t39 = t33 * t36 - t34 * t35;
	t27 = t39 ^ 2;
	t24 = t27 * t28 + 0.1e1;
	t40 = t31 * t43;
	t41 = t39 / t31 * t43;
	t42 = (-t27 * t41 - t39 * t40) / t24 ^ 2;
	t22 = 0.1e1 / t24;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, -0.2e1 * t42 - 0.2e1 * (t22 * t40 - (-t22 * t41 - t28 * t42) * t39) * t39, 0;];
	JaD_rot = t1;
elseif link_index == 4
	%% Symbolic Calculation
	% From jacobiaD_rot_4_floatb_twist_matlab.m
	% OptimizationMode: 2
	% StartTime: 2019-10-09 20:23:31
	% EndTime: 2019-10-09 20:23:31
	% DurationCPUTime: 0.10s
	% Computational Cost: add. (260->9), mult. (248->20), div. (36->4), fcn. (280->4), ass. (0->17)
	t50 = qJ(3) + qJ(4);
	t47 = sin(t50);
	t48 = cos(t50);
	t51 = sin(pkin(6));
	t52 = cos(pkin(6));
	t45 = t47 * t51 + t48 * t52;
	t42 = 0.1e1 / t45 ^ 2;
	t58 = t42 * (qJD(3) + qJD(4));
	t44 = t47 * t52 - t51 * t48;
	t41 = t44 ^ 2;
	t38 = t41 * t42 + 0.1e1;
	t55 = t45 * t58;
	t56 = t44 / t45 * t58;
	t57 = (t41 * t56 + t44 * t55) / t38 ^ 2;
	t36 = 0.1e1 / t38;
	t34 = -0.2e1 * t57 + 0.2e1 * (t36 * t55 + (t36 * t56 - t42 * t57) * t44) * t44;
	t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, t34, t34;];
	JaD_rot = t1;
else
	JaD_rot=NaN(3,4);
end