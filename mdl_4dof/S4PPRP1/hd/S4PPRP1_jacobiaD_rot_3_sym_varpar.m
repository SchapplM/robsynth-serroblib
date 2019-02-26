% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4PPRP1
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
%   pkin=[a2,a3,a4,d3,theta1]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPRP1_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_jacobiaD_rot_3_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_jacobiaD_rot_3_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_jacobiaD_rot_3_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:19:19
% EndTime: 2019-02-26 19:19:19
% DurationCPUTime: 0.03s
% Computational Cost: add. (46->7), mult. (124->20), div. (18->4), fcn. (140->4), ass. (0->15)
t33 = sin(pkin(5));
t34 = cos(pkin(5));
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
JaD_rot  = t1;
