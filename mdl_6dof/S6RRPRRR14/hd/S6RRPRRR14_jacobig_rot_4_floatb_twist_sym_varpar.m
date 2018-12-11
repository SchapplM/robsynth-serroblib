% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:38
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRRR14_jacobig_rot_4_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_jacobig_rot_4_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_jacobig_rot_4_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-10 18:38:21
% EndTime: 2018-12-10 18:38:21
% DurationCPUTime: 0.06s
% Computational Cost: add. (64->27), mult. (69->41), div. (0->0), fcn. (78->19), ass. (0->29)
t210 = sin(pkin(6));
t215 = sin(qJ(1));
t219 = t215 * t210;
t217 = cos(qJ(1));
t218 = t217 * t210;
t216 = cos(qJ(2));
t214 = sin(qJ(2));
t213 = cos(pkin(6));
t212 = cos(pkin(7));
t211 = cos(pkin(8));
t209 = sin(pkin(7));
t208 = sin(pkin(8));
t207 = sin(pkin(14));
t206 = pkin(6) - qJ(2);
t205 = pkin(6) + qJ(2);
t204 = pkin(7) - pkin(14);
t203 = pkin(7) + pkin(14);
t202 = cos(t205);
t201 = sin(t206);
t200 = cos(t206) / 0.2e1;
t199 = sin(t205) / 0.2e1;
t198 = t200 + t202 / 0.2e1;
t197 = t199 - t201 / 0.2e1;
t196 = t199 + t201 / 0.2e1;
t195 = cos(t204) / 0.2e1 + cos(t203) / 0.2e1;
t194 = sin(t203) / 0.2e1 + sin(t204) / 0.2e1;
t193 = -t215 * t198 - t217 * t214;
t192 = t217 * t198 - t215 * t214;
t1 = [0, t219, 0 -(-(-t215 * t197 + t217 * t216) * t207 + t193 * t195 + t194 * t219) * t208 + (-t193 * t209 + t212 * t219) * t211, 0, 0; 0, -t218, 0 -(-(t217 * t197 + t215 * t216) * t207 + t192 * t195 - t194 * t218) * t208 + (-t192 * t209 - t212 * t218) * t211, 0, 0; 1, t213, 0 -(-(t200 - t202 / 0.2e1) * t207 + t196 * t195 + t213 * t194) * t208 + (-t196 * t209 + t213 * t212) * t211, 0, 0;];
Jg_rot  = t1;
