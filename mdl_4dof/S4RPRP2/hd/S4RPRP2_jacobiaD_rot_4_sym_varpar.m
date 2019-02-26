% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
%
% Output:
% JaD_rot [3x4]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S4RPRP2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_jacobiaD_rot_4_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_jacobiaD_rot_4_sym_varpar: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_jacobiaD_rot_4_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:33:01
% EndTime: 2019-02-26 19:33:01
% DurationCPUTime: 0.08s
% Computational Cost: add. (128->13), mult. (322->30), div. (40->4), fcn. (356->4), ass. (0->21)
t56 = sin(qJ(3));
t57 = cos(qJ(1));
t69 = sin(qJ(1));
t70 = cos(qJ(3));
t61 = t56 * t69 + t57 * t70;
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
JaD_rot  = t1;
