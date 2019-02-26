% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PPRRPR1_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_jacobigD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:40:29
% EndTime: 2019-02-26 19:40:29
% DurationCPUTime: 0.04s
% Computational Cost: add. (10->10), mult. (40->28), div. (0->0), fcn. (44->10), ass. (0->14)
t177 = sin(pkin(11));
t183 = cos(pkin(6));
t188 = t177 * t183;
t178 = sin(pkin(7));
t179 = sin(pkin(6));
t187 = t179 * t178;
t181 = cos(pkin(11));
t186 = t181 * t183;
t185 = cos(qJ(3));
t184 = sin(qJ(3));
t182 = cos(pkin(7));
t180 = cos(pkin(12));
t176 = sin(pkin(12));
t1 = [0, 0, 0 ((-t176 * t188 + t181 * t180) * t185 + ((-t181 * t176 - t180 * t188) * t182 + t177 * t187) * t184) * qJD(3), 0, 0; 0, 0, 0 ((t176 * t186 + t177 * t180) * t185 + ((-t177 * t176 + t180 * t186) * t182 - t181 * t187) * t184) * qJD(3), 0, 0; 0, 0, 0 (t178 * t183 * t184 + (t180 * t182 * t184 + t176 * t185) * t179) * qJD(3), 0, 0;];
JgD_rot  = t1;
