% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6PRRPRP1_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:01:11
% EndTime: 2019-02-26 20:01:11
% DurationCPUTime: 0.08s
% Computational Cost: add. (15->10), mult. (24->20), div. (0->0), fcn. (40->8), ass. (0->14)
t103 = sin(pkin(10));
t104 = sin(pkin(6));
t112 = t103 * t104;
t105 = cos(pkin(10));
t111 = t105 * t104;
t106 = cos(pkin(6));
t107 = sin(qJ(2));
t110 = t106 * t107;
t108 = cos(qJ(2));
t109 = t106 * t108;
t102 = qJ(3) + pkin(11);
t101 = cos(t102);
t100 = sin(t102);
t1 = [0, t112, t103 * t109 + t105 * t107, 0 (-t103 * t110 + t105 * t108) * t100 - t101 * t112, 0; 0, -t111, t103 * t107 - t105 * t109, 0 (t103 * t108 + t105 * t110) * t100 + t101 * t111, 0; 0, t106, -t104 * t108, 0, t104 * t107 * t100 - t106 * t101, 0;];
Jg_rot  = t1;
