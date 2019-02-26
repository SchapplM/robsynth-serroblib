% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPPRRR1
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:38
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPPRRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:38:54
% EndTime: 2019-02-26 19:38:54
% DurationCPUTime: 0.10s
% Computational Cost: add. (431->14), mult. (1286->36), div. (18->4), fcn. (1680->14), ass. (0->28)
t139 = sin(pkin(13));
t145 = cos(pkin(13));
t146 = cos(pkin(12));
t140 = sin(pkin(12));
t156 = t140 * cos(pkin(6));
t137 = -t139 * t156 + t146 * t145;
t138 = sin(pkin(14));
t144 = cos(pkin(14));
t136 = -t146 * t139 - t145 * t156;
t142 = sin(pkin(7));
t148 = cos(pkin(7));
t157 = t140 * sin(pkin(6));
t154 = t136 * t148 + t142 * t157;
t133 = t137 * t144 + t154 * t138;
t150 = sin(qJ(4));
t151 = cos(qJ(4));
t155 = (-t136 * t142 + t148 * t157) * sin(pkin(8)) + (-t137 * t138 + t154 * t144) * cos(pkin(8));
t129 = t133 * t150 - t155 * t151;
t130 = t133 * t151 + t155 * t150;
t127 = 0.1e1 / t130 ^ 2;
t164 = qJD(4) * t127;
t126 = t129 ^ 2;
t123 = t126 * t127 + 0.1e1;
t161 = t130 * t164;
t162 = t129 / t130 * t164;
t163 = (t126 * t162 + t129 * t161) / t123 ^ 2;
t121 = 0.1e1 / t123;
t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, 0, -0.2e1 * t163 + 0.2e1 * (t121 * t161 + (t121 * t162 - t127 * t163) * t129) * t129, 0, 0;];
JaD_rot  = t1;
