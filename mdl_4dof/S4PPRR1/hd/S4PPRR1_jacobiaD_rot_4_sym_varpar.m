% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta1]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:21
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4PPRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:21:42
% EndTime: 2019-02-26 19:21:42
% DurationCPUTime: 0.07s
% Computational Cost: add. (260->9), mult. (248->20), div. (36->4), fcn. (280->4), ass. (0->17)
t50 = qJ(3) + qJ(4);
t47 = sin(t50);
t48 = cos(t50);
t51 = sin(pkin(6));
t52 = cos(pkin(6));
t45 = t51 * t47 + t52 * t48;
t42 = 0.1e1 / t45 ^ 2;
t58 = t42 * (qJD(3) + qJD(4));
t44 = t52 * t47 - t51 * t48;
t41 = t44 ^ 2;
t38 = t41 * t42 + 0.1e1;
t55 = t45 * t58;
t56 = t44 / t45 * t58;
t57 = (t41 * t56 + t44 * t55) / t38 ^ 2;
t36 = 0.1e1 / t38;
t34 = -0.2e1 * t57 + 0.2e1 * (t36 * t55 + (t36 * t56 - t42 * t57) * t44) * t44;
t1 = [0, 0, 0, 0; 0, 0, 0, 0; 0, 0, t34, t34;];
JaD_rot  = t1;
