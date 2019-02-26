% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRPR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:35
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->13), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t183 = sin(pkin(6));
t186 = sin(qJ(3));
t196 = t183 * t186;
t187 = sin(qJ(2));
t195 = t183 * t187;
t185 = cos(pkin(6));
t194 = t185 * t187;
t189 = cos(qJ(2));
t193 = t185 * t189;
t192 = qJD(2) * t186;
t182 = sin(pkin(11));
t184 = cos(pkin(11));
t191 = t182 * t189 + t184 * t194;
t190 = -t182 * t194 + t184 * t189;
t188 = cos(qJ(3));
t181 = t183 * t189 * t192 + (t185 * t186 + t188 * t195) * qJD(3);
t180 = (t182 * t196 + t190 * t188) * qJD(3) + (-t182 * t193 - t184 * t187) * t192;
t179 = (-t184 * t196 + t191 * t188) * qJD(3) + (-t182 * t187 + t184 * t193) * t192;
t1 = [0, 0, t190 * qJD(2), t180, 0, -t180; 0, 0, t191 * qJD(2), t179, 0, -t179; 0, 0, qJD(2) * t195, t181, 0, -t181;];
JgD_rot  = t1;
