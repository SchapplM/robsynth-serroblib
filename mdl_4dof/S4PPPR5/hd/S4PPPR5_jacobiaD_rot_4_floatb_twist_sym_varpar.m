% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PPPR5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d4,theta2]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 14:05
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPPR5_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR5_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR5_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR5_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:05:23
% EndTime: 2018-11-14 14:05:23
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
t39 = sin(pkin(5));
t40 = cos(pkin(5));
t44 = sin(qJ(4));
t45 = cos(qJ(4));
t30 = t39 * t45 - t40 * t44;
t26 = 0.1e1 / t30 ^ 2;
t46 = qJD(4) * t26;
t29 = -t39 * t44 - t40 * t45;
t28 = t29 ^ 2;
t23 = t28 * t26 + 0.1e1;
t41 = t30 * t46;
t42 = t29 / t30 * t46;
t43 = (-t28 * t42 - t29 * t41) / t23 ^ 2;
t21 = 0.1e1 / t23;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t43 + 0.2e1 * (-t21 * t41 + (-t21 * t42 - t26 * t43) * t29) * t29;];
JaD_rot  = t1;
