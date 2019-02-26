% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:17
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRRRP4_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP4_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:16:54
% EndTime: 2019-02-26 20:16:54
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t189 = sin(pkin(6));
t192 = sin(qJ(3));
t202 = t189 * t192;
t193 = sin(qJ(2));
t201 = t189 * t193;
t191 = cos(pkin(6));
t200 = t191 * t193;
t195 = cos(qJ(2));
t199 = t191 * t195;
t198 = qJD(2) * t192;
t188 = sin(pkin(11));
t190 = cos(pkin(11));
t197 = t188 * t195 + t190 * t200;
t196 = -t188 * t200 + t190 * t195;
t194 = cos(qJ(3));
t187 = t189 * t195 * t198 + (t191 * t192 + t194 * t201) * qJD(3);
t186 = (t188 * t202 + t196 * t194) * qJD(3) + (-t188 * t199 - t190 * t193) * t198;
t185 = (-t190 * t202 + t197 * t194) * qJD(3) + (-t188 * t193 + t190 * t199) * t198;
t1 = [0, 0, t196 * qJD(2), t186, t186, 0; 0, 0, t197 * qJD(2), t185, t185, 0; 0, 0, qJD(2) * t201, t187, t187, 0;];
JgD_rot  = t1;
