% Rotatorische Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% Jg_rot [3x6]
%   Rotatorische Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg_rot = S6RRPRPR12_jacobig_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_jacobig_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_jacobig_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobig_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:22
% EndTime: 2019-02-26 21:44:22
% DurationCPUTime: 0.03s
% Computational Cost: add. (14->9), mult. (24->18), div. (0->0), fcn. (40->8), ass. (0->16)
t110 = sin(pkin(6));
t113 = sin(qJ(1));
t121 = t113 * t110;
t112 = sin(qJ(2));
t120 = t113 * t112;
t114 = cos(qJ(2));
t119 = t113 * t114;
t115 = cos(qJ(1));
t118 = t115 * t110;
t117 = t115 * t112;
t116 = t115 * t114;
t111 = cos(pkin(6));
t109 = qJ(4) + pkin(11);
t108 = cos(t109);
t107 = sin(t109);
t1 = [0, t121, 0, -t111 * t120 + t116, 0, t107 * t121 - (t111 * t119 + t117) * t108; 0, -t118, 0, t111 * t117 + t119, 0, -t107 * t118 - (-t111 * t116 + t120) * t108; 1, t111, 0, t110 * t112, 0, t110 * t114 * t108 + t111 * t107;];
Jg_rot  = t1;
