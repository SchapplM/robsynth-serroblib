% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR8
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR8_jacobigD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_jacobigD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_jacobigD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:07
% EndTime: 2019-02-26 20:08:07
% DurationCPUTime: 0.06s
% Computational Cost: add. (24->17), mult. (88->43), div. (0->0), fcn. (92->10), ass. (0->24)
t183 = sin(pkin(7));
t184 = sin(pkin(6));
t202 = t184 * t183;
t186 = cos(pkin(7));
t188 = sin(qJ(3));
t201 = t186 * t188;
t187 = cos(pkin(6));
t189 = sin(qJ(2));
t200 = t187 * t189;
t191 = cos(qJ(2));
t199 = t187 * t191;
t198 = t188 * t189;
t190 = cos(qJ(3));
t197 = t190 * t191;
t196 = qJD(2) * t190;
t182 = sin(pkin(12));
t185 = cos(pkin(12));
t195 = -t182 * t189 + t185 * t199;
t194 = t182 * t191 + t185 * t200;
t193 = -t182 * t199 - t185 * t189;
t192 = t182 * t200 - t185 * t191;
t181 = t192 * qJD(2);
t180 = t194 * qJD(2);
t1 = [0, 0, -t181 * t183, 0, t181 * t201 + t193 * t196 + (t192 * t188 + (t182 * t202 + t193 * t186) * t190) * qJD(3), 0; 0, 0, t180 * t183, 0, -t180 * t201 + t195 * t196 + (-t194 * t188 + (-t185 * t202 + t195 * t186) * t190) * qJD(3), 0; 0, 0, qJD(2) * t189 * t202, 0, t187 * t183 * qJD(3) * t190 + ((t186 * t197 - t198) * qJD(3) + (-t186 * t198 + t197) * qJD(2)) * t184, 0;];
JgD_rot  = t1;
