% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR1
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:16:03
% EndTime: 2019-02-26 22:16:03
% DurationCPUTime: 0.02s
% Computational Cost: add. (83->18), mult. (16->4), div. (0->0), fcn. (40->4), ass. (0->10)
t32 = qJ(2) + qJ(3) + pkin(11) + qJ(5);
t30 = sin(t32);
t33 = sin(qJ(1));
t38 = t33 * t30;
t31 = cos(t32);
t37 = t33 * t31;
t34 = cos(qJ(1));
t36 = t34 * t30;
t35 = t34 * t31;
t1 = [-t37, -t36, -t36, 0, -t36, 0; t35, -t38, -t38, 0, -t38, 0; 0, t31, t31, 0, t31, 0; t38, -t35, -t35, 0, -t35, 0; -t36, -t37, -t37, 0, -t37, 0; 0, -t30, -t30, 0, -t30, 0; t34, 0, 0, 0, 0, 0; t33, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
