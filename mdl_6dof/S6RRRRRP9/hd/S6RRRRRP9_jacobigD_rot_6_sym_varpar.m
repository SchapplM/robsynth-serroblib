% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRRRRP9_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP9_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:44:26
% EndTime: 2019-02-26 22:44:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (38->16), mult. (130->40), div. (0->0), fcn. (134->8), ass. (0->24)
t181 = sin(pkin(6));
t186 = cos(qJ(3));
t200 = t181 * t186;
t188 = cos(qJ(1));
t199 = t181 * t188;
t184 = sin(qJ(2));
t185 = sin(qJ(1));
t198 = t184 * t185;
t197 = t184 * t188;
t187 = cos(qJ(2));
t196 = t185 * t187;
t195 = t188 * t187;
t194 = qJD(1) * t181;
t183 = sin(qJ(3));
t193 = qJD(2) * t183;
t182 = cos(pkin(6));
t192 = t182 * t195 - t198;
t191 = t182 * t196 + t197;
t190 = t182 * t197 + t196;
t189 = -t182 * t198 + t195;
t180 = t181 * t187 * t193 + (t182 * t183 + t184 * t200) * qJD(3);
t179 = (-t183 * t199 + t190 * t186) * qJD(3) + t192 * t193 + (t189 * t183 - t185 * t200) * qJD(1);
t178 = (t185 * t181 * t183 + t189 * t186) * qJD(3) - t191 * t193 + (-t190 * t183 - t186 * t199) * qJD(1);
t1 = [0, t188 * t194, t192 * qJD(1) + t189 * qJD(2), t178, t178, 0; 0, t185 * t194, t191 * qJD(1) + t190 * qJD(2), t179, t179, 0; 0, 0, t181 * qJD(2) * t184, t180, t180, 0;];
JgD_rot  = t1;
