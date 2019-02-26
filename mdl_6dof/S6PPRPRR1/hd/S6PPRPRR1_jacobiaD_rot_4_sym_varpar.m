% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRPRR1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.17s
% Computational Cost: add. (337->14), mult. (1021->37), div. (25->4), fcn. (1288->12), ass. (0->26)
t112 = sin(pkin(12));
t117 = cos(pkin(12));
t118 = cos(pkin(11));
t113 = sin(pkin(11));
t129 = t113 * cos(pkin(6));
t108 = -t112 * t129 + t118 * t117;
t111 = sin(pkin(13));
t116 = cos(pkin(13));
t121 = sin(qJ(3));
t122 = cos(qJ(3));
t110 = t121 * t111 - t122 * t116;
t127 = t122 * t111 + t121 * t116;
t138 = sin(pkin(7)) * t113 * sin(pkin(6)) + (-t118 * t112 - t117 * t129) * cos(pkin(7));
t125 = -t108 * t110 + t138 * t127;
t97 = 0.1e1 / t125 ^ 2;
t139 = t125 * t97;
t100 = -t108 * t127 - t110 * t138;
t96 = 0.1e1 / t125;
t94 = t100 * qJD(3);
t133 = t94 * t96 * t97;
t134 = qJD(3) * t139;
t95 = t100 ^ 2;
t92 = t95 * t97 + 0.1e1;
t135 = (-t100 * t134 - t95 * t133) / t92 ^ 2;
t90 = 0.1e1 / t92;
t1 = [0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0; 0, 0, -0.2e1 * t125 * t96 * t135 + 0.2e1 * (-t90 * t134 + (-t90 * t133 - t97 * t135) * t100) * t100 + (t96 - t139) * t90 * t94, 0, 0, 0;];
JaD_rot  = t1;
