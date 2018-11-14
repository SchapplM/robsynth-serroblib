% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:38
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_jacobiaD_rot_4_floatb_twist_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:38:30
% EndTime: 2018-11-14 13:38:30
% DurationCPUTime: 0.04s
% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
t36 = cos(pkin(5));
t37 = cos(qJ(4));
t42 = sin(pkin(5));
t46 = sin(qJ(4));
t31 = t36 * t46 + t42 * t37;
t27 = 0.1e1 / t31 ^ 2;
t47 = qJD(4) * t27;
t30 = t36 * t37 - t42 * t46;
t29 = t30 ^ 2;
t24 = t29 * t27 + 0.1e1;
t43 = t31 * t47;
t44 = t30 / t31 * t47;
t45 = (-t29 * t44 - t30 * t43) / t24 ^ 2;
t22 = 0.1e1 / t24;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t45 + 0.2e1 * (-t22 * t43 + (-t22 * t44 - t27 * t45) * t30) * t30;];
JaD_rot  = t1;
