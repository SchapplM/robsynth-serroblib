% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRP2_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:51:05
% EndTime: 2019-02-26 19:51:05
% DurationCPUTime: 0.06s
% Computational Cost: add. (27->14), mult. (97->37), div. (0->0), fcn. (104->10), ass. (0->19)
t197 = sin(pkin(11));
t200 = cos(pkin(11));
t204 = sin(qJ(2));
t206 = cos(qJ(2));
t208 = t204 * t197 - t206 * t200;
t211 = t208 * qJD(2);
t209 = t197 * t206 + t200 * t204;
t195 = t209 * qJD(2);
t199 = sin(pkin(6));
t203 = sin(qJ(4));
t210 = t199 * t203;
t205 = cos(qJ(4));
t202 = cos(pkin(6));
t201 = cos(pkin(10));
t198 = sin(pkin(10));
t193 = t209 * t202;
t192 = t202 * t211;
t191 = t202 * t195;
t1 = [0, 0, 0, -t198 * t191 - t201 * t211 (t198 * t192 - t201 * t195) * t203 + ((-t198 * t193 - t201 * t208) * t205 + t198 * t210) * qJD(4), 0; 0, 0, 0, t201 * t191 - t198 * t211 (-t201 * t192 - t198 * t195) * t203 + ((t201 * t193 - t198 * t208) * t205 - t201 * t210) * qJD(4), 0; 0, 0, 0, t199 * t195, t202 * qJD(4) * t203 + (t209 * qJD(4) * t205 - t203 * t211) * t199, 0;];
JgD_rot  = t1;
