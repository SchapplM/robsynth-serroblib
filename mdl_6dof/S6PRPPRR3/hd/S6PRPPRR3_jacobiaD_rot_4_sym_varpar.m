% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPPRR3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:45:50
% EndTime: 2019-02-26 19:45:50
% DurationCPUTime: 0.10s
% Computational Cost: add. (153->13), mult. (451->29), div. (25->4), fcn. (552->7), ass. (0->21)
t73 = cos(pkin(10));
t75 = sin(qJ(2));
t76 = cos(qJ(2));
t82 = sin(pkin(10)) * cos(pkin(6));
t67 = t73 * t75 + t76 * t82;
t68 = t73 * t76 - t75 * t82;
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t62 = t67 * t70 + t68 * t72;
t57 = 0.1e1 / t62 ^ 2;
t79 = -t67 * t72 + t68 * t70;
t89 = t57 * t79 ^ 2;
t65 = t67 * qJD(2);
t66 = t68 * qJD(2);
t54 = -t65 * t72 + t66 * t70;
t56 = 0.1e1 / t62;
t86 = t54 * t56;
t88 = t86 * t89;
t85 = t79 * (-t65 * t70 - t66 * t72);
t52 = 0.1e1 + t89;
t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0.2e1 * (t56 * t62 + t89) / t52 ^ 2 * (t57 * t85 - t88) + (-t86 + 0.2e1 * t88 + (t54 * t62 - 0.2e1 * t85) * t57) / t52, 0, 0, 0, 0;];
JaD_rot  = t1;
