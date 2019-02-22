% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 10:17
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RPPRRR1_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 10:17:35
% EndTime: 2019-02-22 10:17:35
% DurationCPUTime: 0.02s
% Computational Cost: add. (58->14), mult. (12->4), div. (0->0), fcn. (30->4), ass. (0->11)
t27 = pkin(11) + qJ(4) + qJ(5);
t23 = sin(t27);
t28 = qJ(1) + pkin(10);
t25 = sin(t28);
t32 = t25 * t23;
t24 = cos(t27);
t31 = t25 * t24;
t26 = cos(t28);
t30 = t26 * t23;
t29 = t26 * t24;
t1 = [-t31, 0, 0, -t30, -t30, 0; t29, 0, 0, -t32, -t32, 0; 0, 0, 0, t24, t24, 0; t32, 0, 0, -t29, -t29, 0; -t30, 0, 0, -t31, -t31, 0; 0, 0, 0, -t23, -t23, 0; t26, 0, 0, 0, 0, 0; t25, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0;];
JR_rot  = t1;
