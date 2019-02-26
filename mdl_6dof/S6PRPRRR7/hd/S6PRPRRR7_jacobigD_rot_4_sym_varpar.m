% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR7_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobigD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobigD_rot_4_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:30
% EndTime: 2019-02-26 19:57:30
% DurationCPUTime: 0.04s
% Computational Cost: add. (12->10), mult. (54->29), div. (0->0), fcn. (54->12), ass. (0->15)
t205 = sin(pkin(7)) * cos(pkin(8));
t204 = cos(pkin(7)) * cos(pkin(14));
t198 = cos(pkin(6));
t199 = sin(qJ(2));
t203 = t198 * t199;
t200 = cos(qJ(2));
t202 = t198 * t200;
t189 = sin(pkin(14));
t201 = qJD(2) * t189;
t195 = cos(pkin(13));
t191 = sin(pkin(8));
t190 = sin(pkin(13));
t188 = (t190 * t203 - t195 * t200) * qJD(2);
t187 = (-t190 * t200 - t195 * t203) * qJD(2);
t1 = [0, 0, 0 -(t188 * t204 - (-t190 * t202 - t195 * t199) * t201) * t191 - t188 * t205, 0, 0; 0, 0, 0 -(t187 * t204 - (-t190 * t199 + t195 * t202) * t201) * t191 - t187 * t205, 0, 0; 0, 0, 0 (-(-t189 * t200 - t199 * t204) * t191 + t199 * t205) * sin(pkin(6)) * qJD(2), 0, 0;];
JgD_rot  = t1;
